import os
from os.path import join, splitext, exists, isdir, isfile
import subprocess
import logging

from fs.opener import opener
# from fs.zipfs import ZipFS


logger = logging.getLogger()


def run_fastqc(fastq_paths,
               output_directory=None,
               fastqc_bin=None,
               threads=2,
               extra_options=''):

    if not fastq_paths:
        logger.warning('FastQC - called with no FASTQ file paths provided, '
                       'skipping.')
        return None

    tmp_dir = None
    if not output_directory:
        try:
            tmp_dir = create_tmp_dir()
            output_directory = tmp_dir
        except OSError:
            logger.error('FastQC - failed to create temp directory.')
            return None

    if not fastqc_bin:
        # We will assume it's on path with the default name
        fastqc_bin = 'fastqc'

    cmd = 'nice %s %s --noextract --threads %d --outdir %s %s' % \
          (fastqc_bin,
           extra_options,
           threads,
           output_directory,
           ' '.join(fastq_paths))

    logger.info('Running FastQC on: %s', ', '.join(fastq_paths))

    cmd_out = None
    try:
        # Unfortunately FastQC doesn't always return sensible
        # exit codes on failure, so we parse the output
        cmd_out = subprocess.check_output(cmd,
                                          shell=True,
                                          stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        logger.error('FastQC stdout: %s', cmd_out)
        return None

    success = False
    if 'Analysis complete' in cmd_out.splitlines()[-1]:
        success = True
    else:
        logger.error('FastQC stdout: %s', cmd_out)

    if not success:
        if tmp_dir is not None:
            shutil.rmtree(output_directory)
            return None

    return output_directory


def file_from_zip(zip_file_path, filename, mode='r'):
    # eg rootpath =
    # FastQC.out/15-02380-CE11-T13-L1_AACCAG_L001_R1_001_fastqc.zip
    # ie   fastq_filename + '_fastqc.zip'

    # fs.opener.opener magic detects a filesystem type if we give it
    # a URI-like string, eg zip://foo/bla/file.zip
    # (eg https://commons.apache.org/proper/commons-vfs/filesystems.html)
    # So we munge the path to a zip file to become a compatible URI
    # (an http, ftp, webdav url could also be used)
    if splitext(zip_file_path)[1] == '.zip':
        zip_file_path = 'zip://' + zip_file_path

    with opener.parse(zip_file_path)[0] as vfs:
        for fn in vfs.walkfiles():
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
