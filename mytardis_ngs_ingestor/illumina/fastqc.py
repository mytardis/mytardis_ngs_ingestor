import os
from os.path import join, splitext, exists, isdir, isfile
import shutil
from pathlib2 import Path
import subprocess
import psutil
import logging

try:
    from shutil import which
except ImportError:
    from backports.shutil_which import which

from fs import open_fs
# from fs.zipfs import ZipFS

from mytardis_ngs_ingestor.utils import tmp_dirs

_fastqc_jvm_Xmx_memory = 250  # MB


def _fastq_is_empty(fqfile):
    """
    Return True if a FASTQ files contains no reads.

    :param fqfile: The path to the FASTQ file (filename).
    :type fqfile: str
    :return: True if the file contain no reads
    :rtype: bool
    """
    cmd = 'gunzip -c %s | head | wc -c' % fqfile
    logging.debug('Command: %s', cmd)
    try:
        cmd_out = subprocess.check_output(cmd,
                                          shell=True,
                                          stderr=subprocess.STDOUT)
        if cmd_out.strip() == '0':
            logging.info('File %s contains no reads - skipping FastQC on this file.', fqfile)
            return True

    except subprocess.CalledProcessError:
        logging.warn('Failed attempting to check for empty FASTQ file: %s', cmd)

    return False


def run_fastqc(fastq_paths,
               output_directory=None,
               fastqc_bin=None,
               threads=None,
               extra_options=''):
    cmd_out = None

    if not fastqc_bin:
        # We will assume it's on path with the default name
        fastqc_bin = 'fastqc'

    if not which(fastqc_bin):
        logging.error("FastQC - executable %s doesn't exist ?", fastqc_bin)
        return None, cmd_out

    if threads is None:
        # We set threads to the number of cpus, or less if the RAM would be
        # exhausted
        cpus = psutil.cpu_count()
        threads = int(cpus)

        memory = psutil.virtual_memory().total / (1024.0 * 1024.0)  # MB
        max_threads_for_mem = memory / _fastqc_jvm_Xmx_memory
        if max_threads_for_mem < cpus:
            threads = int(max_threads_for_mem) - 1

    if not fastq_paths:
        logging.warning('FastQC - called with no FASTQ file paths provided, '
                        'skipping.')
        return None, cmd_out

    for f in fastq_paths:
        if _fastq_is_empty(f):
            fastq_paths.remove(f)

    tmp_dir = None
    if not output_directory:
        try:
            tmp_dir = tmp_dirs.create_tmp_dir()
            output_directory = tmp_dir
        except OSError:
            logging.error('FastQC - failed to create temp directory.')
            return None, cmd_out

    cmd = 'nice %s %s --noextract --threads %d --outdir %s %s' % \
          (fastqc_bin,
           extra_options,
           threads,
           output_directory,
           ' '.join(fastq_paths))

    logging.info('Running FastQC on: %s', ', '.join(fastq_paths))
    logging.info('Command: %s', cmd)

    try:
        # Unfortunately FastQC doesn't always return sensible
        # exit codes on failure, so we parse the output
        cmd_out = subprocess.check_output(cmd,
                                          shell=True,
                                          stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        logging.error('FastQC stdout: %s', cmd_out)
        return None, cmd_out

    success = False
    if 'Analysis complete' in cmd_out.splitlines()[-1]:
        success = True
    else:
        logging.error('FastQC stdout: %s', cmd_out)

    if not success:
        if tmp_dir is not None:
            shutil.rmtree(output_directory)
            return None, cmd_out

    return output_directory, cmd_out


def get_fastqc_output_directory(proj_path):
    return join(proj_path, 'FastQC.out')


def run_fastqc_on_project(fastq_files,
                          proj_path,
                          output_directory=None,
                          fastqc_bin=None,
                          threads=2,
                          clobber=False):
    cmd_out = None

    if output_directory is None:
        output_directory = get_fastqc_output_directory(proj_path)
        if exists(output_directory):
            if clobber:
                logging.warning("Removing old FastQC output directory: %s",
                                output_directory)
                shutil.rmtree(output_directory)
            else:
                logging.warning("FastQC - output directory already exists: %s",
                                output_directory)
        else:
            try:
                # os.mkdir(output_directory)
                Path(str(output_directory)).mkdir(parents=True, exist_ok=True)
            except OSError:
                logging.error("FastQC - couldn't create output directory: %s",
                              output_directory)
                return None, cmd_out

    if not exists(output_directory) or \
            not isdir(output_directory):
        logging.error("FastQC - output path %s isn't a directory "
                      "or doesn't exist.",
                      output_directory)
        return None, cmd_out

    fqc_output_directory, cmd_out = run_fastqc(
        fastq_files,
        output_directory=output_directory,
        fastqc_bin=fastqc_bin,
        threads=threads)

    return fqc_output_directory, cmd_out


def file_from_zip(zip_file_path, filename, mode='r'):
    # eg rootpath =
    # FastQC.out/15-02380-CE11-T13-L1_AACCAG_L001_R1_001_fastqc.zip
    # ie   fastq_filename + '_fastqc.zip'

    # fs.open_fs magic detects a filesystem type if we give it
    # a URI-like string, eg zip://foo/bla/file.zip
    # (eg https://commons.apache.org/proper/commons-vfs/filesystems.html)
    # So we munge the path to a zip file to become a compatible URI
    # (an http, ftp, webdav url could also be used)
    if splitext(zip_file_path)[1] == '.zip':
        zip_file_path = 'zip://' + zip_file_path

    with open_fs(zip_file_path) as vfs:
        for fn in vfs.walk.files():
            if os.path.basename(fn) == filename:
                return vfs.open(fn, mode)


def parse_file_from_zip(zip_file_path, filename, parser):
    return parser(file_from_zip(zip_file_path, filename))


def parse_summary_txt(zip_file_path):
    """
    Extract the overall PASS/WARN/FAIL summary.txt table from FastQC results
    contained in a zip file.
    """

    def _parse(fh):
        summary = [tuple(line.strip().split('\t')) for line in fh]
        return summary

    qc_summary = parse_file_from_zip(zip_file_path,
                                     'summary.txt',
                                     _parse)
    return qc_summary


def parse_data_txt(zip_file_path):
    """
    Extract the tables in fastqc_data.txt from FastQC results contained
    in a zip file.

    Returns a dictionary representing the tables in the file,
    keyed by section headers in the form:

     {u'Basic Statistics':
      {'column_labels': (u'Measure', u'Value'),
      'rows': [
       (u'Filename', u'14-06207-SOX-5_AGTGAG_L001_R1_001.fastq.gz'),
       (u'File type', u'Conventional base calls'),
       (u'Encoding', u'Sanger / Illumina 1.9'), (u'Total Sequences', u'157171'),
       (u'Filtered Sequences', u'0'), (u'Sequence length', u'51'),
       (u'%GC', u'40')
      ],
      'qc_result': u'pass'},

      u'fastqc_version': '0.12',
     }

     The FastQC version is stored under data['fastqc_version'].

    :type zip_file_path: str
    :return: dict
    """

    def _parse(fh):
        data = {'fastqc_version': fh.readline().split('\t')[1].strip()}
        section = None
        for l in fh:
            line = l.strip()
            if line[0:2] == '>>':
                if line == '>>END_MODULE':
                    # end section
                    continue
                else:
                    # start new section
                    section, qc_result = line[2:].split('\t')
                    data[section] = {'rows': [],
                                     'qc_result': qc_result}
                    continue
            if line[0] == '#':
                column_labels = tuple(line[1:].split('\t'))
                data[section]['column_labels'] = column_labels
                continue
            else:
                data_row = tuple(line.split('\t'))
                data[section]['rows'].append(data_row)

        return data

    qc_data = parse_file_from_zip(zip_file_path,
                                  'fastqc_data.txt',
                                  _parse)

    return qc_data


def extract_basic_stats(fastqc_data, sample_id):
    basic_stats = fastqc_data[sample_id]['Basic Statistics']
    stats = {}
    for row in basic_stats['rows']:
        k = row[0]
        v = row[1]
        if 'Total Sequences' in k:
            stats['number_of_reads'] = int(v)
        # No longer used since apparently this is always zero ?
        # if 'Sequences flagged as poor quality' in k:
        #     stats['number_of_poor_quality_reads'] = int(v)
        if 'Sequence length' in k:
            # this can be a single number (eg 51) or a range (eg 35-51)
            # we take the upper value in the range
            stats['read_length'] = int(v.split('-')[-1])
        if '%GC' in k:
            stats['percent_gc'] = float(v)

    return stats
