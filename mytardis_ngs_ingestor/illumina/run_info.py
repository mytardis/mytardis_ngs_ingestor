from __future__ import absolute_import, division, print_function
from builtins import (bytes, str, open, super, range,
                      zip, round, input, int, pow, object)
import logging
import os
from os.path import join, splitext, exists, isdir, isfile
from pathlib2 import Path
from fs import open_fs
import subprocess
import re
import csv
import gzip
from collections import OrderedDict
from dateutil import parser as dateparser
import xmltodict

logger = logging.getLogger()


def parse_samplesheet(file_path, standardize_keys=True):

    # Old plain CSV format, IEM v3:
    # FCID,Lane,SampleID,SampleRef,Index,Description,Control,Recipe,
    # Operator,SampleProject
    #
    # last line: #_IEMVERSION_3_TruSeq LT,,,,,,,,,

    # Newer INI-style CSV format, IEM v4:
    # [Header],,,,,,,,
    # IEMFileVersion,4,,,,,,,
    # .. >snip< ..
    # Assay,TruSeq LT,,,,,,,
    # [Reads],,,,,,,,
    # .. >snip< ..
    # [Settings],,,,,,,,
    # .. >snip< ..
    # [Data],,,,,,,,
    # Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,Sample_Project,Description
    # .. >snip< ..
    #

    def _format_keys(reader):
        if standardize_keys:
            formatted = []
            # remove any underscores to make names consistent between
            # old and new style samplesheets (eg Sample_ID -> SampleID)
            for r in [row for row in reader]:
                r = {k.replace('_', '').lower(): r[k] for k in r.keys()}
                formatted.append(r)
        else:
            formatted = [row for row in reader]

        return formatted

    with open(file_path, "rU") as f:
        lines = f.readlines()

    chemistry = None

    # IEM v4 INI-style CSV
    if '[Header]' in lines[0]:
        section = None
        for i, l in enumerate(lines):
            if l[0] == '[':
                section = l[1:].split(']')[0]
            if section == 'Header' and l.startswith('Assay,'):
                chemistry = l.split(',')[1].strip()
            if section == 'Data':
                data_index = i
                break

        reader = csv.DictReader(lines[data_index+1:])

        samples = _format_keys(reader)

    else:  # Plain CSV (IEM v3 ?)
        reader = csv.DictReader(lines)
        samples = _format_keys(reader)
        lastlinebit = samples[-1:][0].get('FCID', None) or \
                      samples[-1:][0].get('fcid', None)
        if lastlinebit is not None:
            chemistry = lastlinebit.split('_')[-1].strip()
        del samples[-1:]

    return samples, chemistry


def filter_samplesheet_by_project(file_path, proj_id,
                                  project_column_label='SampleProject',
                                  output_ini_headers=False):

    """
    Windows \r\n

    :param file_path:
    :type file_path:
    :param proj_id:
    :type proj_id:
    :param project_column_label:
    :type project_column_label:
    :return:
    :rtype:
    """
    # FCID,Lane,SampleID,SampleRef,Index,Description,Control,Recipe,
    # Operator,SampleProject
    #
    # last line: #_IEMVERSION_3_TruSeq LT,,,,,,,,,
    outlines = []
    with open(file_path, "rU") as f:
        header_content = []
        header = f.readline().strip()
        header_content.append(header + '\r\n')

        # skip any INI-style headers down to the CSV sample list in the [Data]
        if '[Header]' in header:
            while '[Data]' not in header:
                header = f.readline().strip()
                header_content.append(header + '\r\n')
            header = f.readline().strip()
            header_content.append(header + '\r\n')

        if output_ini_headers:
            outlines.extend(header_content[:-1])

        s = header.split(',')
        # old samplesheet formats have no underscores in column labels,
        # newer (IEMv4) ones do. by removing underscores here, we can find
        # 'SampleProject' and 'Sample_Project', whichever exists
        s_no_underscores = [c.lower().replace('_', '') for c in s]
        project_column_index = s_no_underscores.index(project_column_label.lower())
        outlines.append(header + '\r\n')

        # find the set of unique project ids in the sample sheet
        proj_ids_in_sheet = set()
        for l in f:
            s = l.strip().split(',')
            proj_ids_in_sheet.add(s[project_column_index].strip())

        # We want to output all lines if there is only a single project.
        # The on-instrument FASTQ generation for MiSeq seems to not produce
        # Project directories, even if there is a Sample_Project specified in
        # the samplesheet. So the 'proj_id' here might be empty, but we still
        # want all the SampleSheet.csv lines of the one and only project in the
        # SampleSheet in this case.
        one_project_only = (len(proj_ids_in_sheet) == 1)

        for l in f:
            s = l.strip().split(',')
            if (one_project_only or
                    s[project_column_index].strip() == proj_id.strip() or
                    l[0] == '#'):
                outlines.append(l.strip() + '\r\n')
    return outlines


def samplesheet_to_dict(samplesheet_rows, key='SampleID'):
    by_sample_id = {}
    key = key.lower()
    for sample in samplesheet_rows:
        sample_id = sample[key]
        by_sample_id[sample_id] = sample.copy()
        del by_sample_id[sample_id][key]

    return by_sample_id


def get_project_ids_from_samplesheet(samplesheet, include_no_index_name=None):
    # Get unique project names
    projects = list(set([s['SampleProject'] for s in samplesheet]))

    # Treat the 'Undetermined_indices' directory as a (special) project
    if include_no_index_name is not None:
        projects.append(include_no_index_name)

    return projects


def get_number_of_reads_fastq(filepath):
    """
    Count the number of reads in a (gzipped) FASTQ file.
    Assumes fours lines per read.

    :type filepath: str
    :rtype: int
    """
    num = subprocess.check_output("zcat %s | echo $((`wc -l`/4))" % filepath,
                                  shell=True)
    return int(num.strip())


def get_read_length_fastq(filepath):
    """
    Return the length of the first read in a (gzipped) FASTQ file.

    :param filepath: Path to the (gzipped) FASTQ file
    :type filepath: str
    :return: Length of the first read in the FASTQ file
    :rtype: int
    """
    num = subprocess.check_output("zcat < %s | "
                                  "head -n 2  | "
                                  "tail -n 1  | "
                                  "wc -m" % filepath,
                                  shell=True)
    return int(num.strip()) - 1


# Copypasta from: https://goo.gl/KpWo1w
# def unique(seq):
#     seen = set()
#     seen_add = seen.add
#     return [x for x in seq if not (x in seen or seen_add(x))]


def rta_complete_parser(run_path):
    """
    Parses RTAComplete.txt files in completed Illumina runs.
    Returns the Illumina RTA (Realtime Analysis) version and the
    date & time the run finished transferring from the instrument.

    Copes with file generated by both RTA 1.x and RTA 2.x.

    :type run_path: str
    :rtype: (datetime.DateTime, str)
    """
    with open(join(run_path, "RTAComplete.txt"), 'r') as f:
        line = f.readline()
        if line.startswith('RTA'):
            # RTA 2.7.3 completed on 3/25/2016 3:31:22 AM
            s = line.split(' completed on ')
            version = s[0]
            day, time = s[1].split(' ', 1)
        else:
            # 6/11/2014,20:00:49.935,Illumina RTA 1.17.20
            day, time, version = line.split(',')

    end_time = dateparser.parse("%s %s" % (day, time))
    return end_time, version


def runinfo_parser(run_path):
    """

    Matches some or all of the fields defined in schema:
    http://www.tardis.edu.au/schemas/sequencing/run/illumina

    :type run_path: str
    :rtype: dict
    """
    with open(join(run_path, "RunInfo.xml"), 'r') as f:
        runinfo = xmltodict.parse(f.read())['RunInfo']['Run']

    info = {u'run_id': runinfo['@Id'],
            u'run_number': runinfo['@Number'],
            u'flowcell_id': runinfo['Flowcell'],
            u'instrument_id': runinfo['Instrument']}

    reads = runinfo['Reads']['Read']

    cycle_list = []
    # index_reads = []
    for read in reads:
        if read['@IsIndexedRead'] == 'Y':
            # index_reads.append(read['@Number'])
            # we wrap the index reads in brackets
            cycle_list.append("(%s)" % read['@NumCycles'])
        else:
            cycle_list.append(read['@NumCycles'])

    info['read_cycles'] = ', '.join(cycle_list)
    # info['index_reads'] = ', '.join(index_reads)

    # Currently not capturing this metadata
    # runinfo['RunInfo']['Run']['FlowcellLayout']['@LaneCount']
    # runinfo['RunInfo']['Run']['FlowcellLayout']['@SurfaceCount']
    # runinfo['RunInfo']['Run']['FlowcellLayout']['@SwathCount']
    # runinfo['RunInfo']['Run']['FlowcellLayout']['@TileCount']

    return info


def illumina_config_parser(run_path):
    """
    Extacts data from an Illumina run Config/*_Effective.cfg file.

    Returns a dictionary with key/values, where keys have
    the [section] from the Windows INI-style file are prepended,
    separated by colon. eg

    {"section_name:variable_name": "value"}

    :type run_path: str
    :rtype: dict
    """

    # we find the approriate config file
    config_filename = None
    for filename in os.listdir(join(run_path, 'Config')):
        if "Effective.cfg" in filename:
            config_filename = filename
    if not config_filename:
        logger.error("Cannot find Config/*Effective.cfg file")
        return

    # we don't use ConfigParser since for whatever reason it can't handle
    # these files, despite looking like a mostly sane Windows INI-style
    allinfo = {}
    section = None
    with open(join(run_path, 'Config', config_filename), 'r') as f:
        for l in f:
            if l[0] == ';':
                continue
            if l[0] == '[':
                section = l[1:].split(']')[0]
            if '=' in l:
                s = l.split('=')
                k = s[0].strip()
                v = s[1]
                if ';' in l:
                    v = s[1].split(';')[0]
                v = v.strip()
                allinfo['%s:%s' % (section, k)] = v

    # select just the keys of interest to return
    # info = {}
    # if 'system:instrumenttype' in allinfo:
    #    info['instrument_model'] = allinfo['system:instrumenttype']

    return allinfo


# TODO: Better demultiplexer version & commandline detection
# since we don't have a really reliable way of guessing how the demultiplexing
# was done (beyond detecting DemultiplexConfig.xml for bcl2fastq 1.8.4),
# we should probably require the processing pipeline to output some
# an optional metadata file containing the version and commandline used
# (and possibly other stuff)
# Alternative - allow config / commandline to specify pattern/path to
# bcl2fastq (or other demultiplexer) stdout log - extract version and
# commandline from there. This should work for bcl2fastq 2.17 but not
# 1.8.4
def get_demultiplexer_info(demultiplexed_output_path):
    """
    Determine which demultiplexing program (eg bcl2fastq) and commandline
    options that were used to partition reads from the run into individual
    samples (based on index).

    eg. {'version': 'bcl2fastq 1.8.3,
         'commandline_options':
         '--input-dir ./Data/Intensities/BaseCalls ' +
         '--output-dir ./130613_DMO177_0029_AH0EPTADXX.pc ' +
         '--sample-sheet ./SampleSheet.csv --no-eamss'}

    :param demultiplexed_output_path: bcl2fastq (or similar) output path
    :type demultiplexed_output_path: str
    :return: Name and version number of program used for demultiplexing reads.
    :rtype: dict
    """

    version_info = {'version': '',
                    'version_number': '',
                    'commandline_options': ''}

    # Parse DemultiplexConfig.xml to extract the bcl2fastq version
    # This works for bcl2fastq v1.8.4, but bcl2fastq2 v2.x doesn't seem
    # to generate this file
    demulti_config_path = join(demultiplexed_output_path,
                               "DemultiplexConfig.xml")
    if exists(demulti_config_path):
        with open(demulti_config_path, 'r') as f:
            xml = xmltodict.parse(f.read())
            version = xml['DemultiplexConfig']['Software']['@Version']
            version = ' '.join(version.split('-'))
            cmdline = xml['DemultiplexConfig']['Software']['@CmdAndArgs']
            cmdline = cmdline.split(' ', 1)[1]
            version_info = {'version': version,
                            'version_number': version.split()[1].lstrip('v'),
                            'commandline_options': cmdline}
    else:
        # if we can't find DemultiplexConfig.xml, assume the locally installed
        # bcl2fastq2 (v2.x) version was used

        # TODO: This assumption is just misleading for things like MiSeq runs
        #       the might have been demultiplexed on the instrument and
        #       copied over. We probably shouldn't do this, just leave it blank
        #       until we find a better way to determine the demultiplexer
        #       (which might amount to having it specified on the commandline or
        #        in an extra metadata file)
        try:
            out = subprocess.check_output("/usr/local/bin/bcl2fastq --version",
                                          stderr=subprocess.STDOUT,
                                          shell=True).splitlines()
            if len(out) >= 2 and 'bcl2fastq' in out[1]:
                version = out[1].strip()
                version_info = {
                    'version': version,
                    'version_number': version.split()[1].lstrip('v'),
                    'commandline_options': None}
        except subprocess.CalledProcessError:
            pass

    # if we can't determine the bcl2fastq (or other demultiplexer) based
    # on config & log files, or the locally installed version, try to guess if
    # bcl2fastq 1.x or 2.x was used based on 'Project_' prefixes
    if not version_info.get('version'):
        if any([proj_dir.startswith('Project_')
               for proj_dir in os.listdir(demultiplexed_output_path)]):
            version_info['version'] = 'bcl2fastq 1.0unknown'
            version_info['version_number'] = '1.0unknown'
        else:
            version_info['version'] = 'bcl2fastq 2.0unknown'
            version_info['version_number'] = '2.0unknown'

    return version_info
    """
    # get the version of locally installed tagdust

    out = subprocess.check_output("tagdust --version", shell=True)
    if len(out) >= 1 and 'Tagdust' in out[1]:
        version = out[1].strip()
        return {'version': version,
                'commandline_options': None}
    """


def get_run_id_from_path(run_path):
    return os.path.basename(run_path.strip(os.path.sep))


# TODO: We could actually make patterns like these one a config option
# (either raw Python regex with named groups, or write a translator to
#  simplify the syntax so we can write {sample_id}_bla_{index} in the config
#  and have it converted to a regex with named groups internally)
# https://docs.python.org/2/howto/regex.html#non-capturing-and-named-groups
#
# (All these problems of differences between v2.x and v1.8.4 and CSV vs.
#  IEMv4 SampleSheets are begging for a SampleSheet data container
#  object to abstract out differences in sample sheets, and an IlluminaRun
#  object with a list of DemultiplexedProject objects to abstract out
#  differences in directory structure)
def parse_sample_info_from_filename(filepath, suffix='.fastq.gz'):
    """
    Given a filename in one of the common formats generated by
    bcl2fastq (or FastQC), return metadata about the file, eg lane, read number,
    read type, sample name, index etc.

    :param filepath: A filename (can include path)
    :type filepath: basestring
    :param suffix: A suffix to trim prior to extracting values.
    :type suffix: basestring
    :return: A dictionary, eg {'lane': 1, 'read': 2, ...}
    :rtype: dict
    """
    filename = os.path.basename(filepath)

    def change_dict_types(d, keys, map_fn, ignore_missing_keys=True):
        for k in keys:
            if k in d:
                d[k] = map_fn(d[k])
            elif not ignore_missing_keys:
                raise KeyError("%s not found" % k)
        return d

    # Components of regexes to match common fastq.gz filenames output
    # by Illumina software. It's easier and more extensible to combine
    # these than attempt to create and maintain one big regex
    sample_name_re = r'(?P<sample_name>.*)'
    undetermined_sample_name_re = r'(?P<sample_name>.*_Undetermined)'
    index_re = r'(?P<index>[ATGC-]{6,33})'
    # dual_index_re = r'(?P<index>[ATGC]{6,12})-?(?P<index2>[ATGC]{6,12})?'
    lane_re = r'L0{0,3}(?P<lane>\d+)'
    # Usually R1 or R2, can be I1 etc if index reads are output seperately
    read_re = r'(?P<read_type>[RI])(?P<read>\d)'
    # index_read_re = r'I(?P<read>\d)'
    set_number_re = r'(?P<set_number>\d+)'
    sample_number_re = r'S(?P<sample_number>\d+)'
    extension_re = r'%s$' % suffix

    filename_regexes = [
        # Undetermined indices files like this:
        # lane1_Undetermined_L001_R1_001.fastq.gz
        r'_'.join([undetermined_sample_name_re, lane_re, read_re, set_number_re]),

        # bcl2fastq 2.x style filenames:
        # {sample_name}_{sample_number}_L00{lane}_R{read}_001.fastq.gz
        r'_'.join([sample_name_re, sample_number_re, lane_re, read_re,
                   set_number_re]),

        # bcl2fastq 1.8.4 style filenames:
        # {sample_name}_{index}_L00{lane}_R{read}_001.fastq.gz
        r'_'.join([sample_name_re, index_re, lane_re, read_re, set_number_re]),
    ]

    filename_regexes = [r + extension_re for r in filename_regexes]

    # path_regexes = [
    #     r'(?P<project_id>.*)/Sample_(?P<sample_id>.*)/'
    #     r'(?P<project_id>.*)/(?P<sample_id>.*)/',
    #     r'',
    # ]

    # combined_re = r'(%s)' % '|'.join(filename_regexes)
    # combined_re = r'(%s)(%s)' % ('|'.join(path_regexes),
    #                              '|'.join(filename_regexes))

    for regex in filename_regexes:
        m = re.match(regex, filename)
        if m is not None:
            break

    if m is not None:
        d = m.groupdict()
        d = change_dict_types(d, ['sample_number', 'lane', 'read',
                                  'set_number'], int)
        return d

    return None


def get_sample_project_mapping(basepath,
                               samplesheet=None,
                               suffix='.fastq.gz',
                               absolute_paths=False,
                               catch_undetermined=True):
    """
    Given a path containing fastq.gz files, possibily nested in Project/Sample
    directories, return a data structure mapping fastq-samples to projects.

    TODO: The SampleSheet.csv may be used as a hint but is not required.

    :param basepath: Path to directory tree of fastq.gz files - eg, bcl2fastq
                     output directory
    :type basepath: str
    :param catch_undetermined: Map any FASTQ files with 'Undetermined' in the
                               name to a projects 'Undetermined_indices'
    :type catch_undetermined: bool
    :return: Dictionary lists, {project_id : [relative fastq.gz paths]}
    :rtype: OrderedDict
    """

    if not os.path.exists(basepath):
        raise IOError('Path %s does not exist !' % basepath)

    fq_files = []
    with open_fs(basepath) as vfs:
        for fn in vfs.walk.files():
            if fn.endswith(suffix):
                fq_files.append(fn.lstrip('/').lstrip('\\'))

    fq_files = sorted(fq_files)

    project_mapping = OrderedDict()
    for fqpath in fq_files:
        project = ''
        fqfile = fqpath
        parts = Path(fqpath).parts
        if len(parts) == 3:
            project, sample_id, fqfile = map(str, parts)
        if len(parts) == 2:
            project, fqfile = map(str, parts)
        if len(parts) == 1:
            fqfile = str(parts[0])

        if catch_undetermined and 'Undetermined' in fqfile:
            project = u'Undetermined_indices'

        # TODO: we currently don't deal with Project_ prefixes, really
        #       the project ID doesn't include Project_. If we strip
        #       this here, maybe we need to include the project directory
        #       in the fastq paths so we can know the path and project id
        #       - will require fixes to downstream code that
        #       does join(bcl2fastq_output_dir, project_id, fastq_file)

        # TODO: also incorporate sample_id in this datastructure
        if project not in project_mapping:
            project_mapping[project] = []
        if absolute_paths:
            fqpath = join(basepath, fqpath)
        project_mapping[project].append(fqpath)

    # TODO: Use the SampleSheet.csv to validate or hint
    # TODO: Also assign sample_id, sample_name, lane, read, etc
    #       we could use parse_sample_info_from_filename for this,
    #       and/or use the FASTQ header(s)

    return project_mapping


def undetermined_reads_in_root(basepath):
    """
    Returns False if the Undetermined_indicies fastq.gz files are in their own
    project directory (ie, old bcl2fastq 1.8.4 style), or True if they are
    bare in the root of the bcl2fastq output diretory (eg bcl2fastq 2.x)

    :param basepath:
    :type basepath:
    :return:
    :rtype:
    """
    for f in os.listdir(basepath):
        if 'Undetermined_' not in f:
            continue

        f_path = join(basepath, f)

        if isfile(f_path):
            return True
        if isdir(f_path):
            return False

    raise Exception("No Undetermined_indices files or directories found.")


def get_sample_id_from_fastq_filename(filepath):
    """
    Takes:
    15-02380-CE11-T13-L1_AACCAG_L001_R1_001.fastq.gz

    Returns a sample ID that should be unique within a run,
    consisting of the sample name, index, lane and read pair:
    15-02380-CE11-T13-L1_AACCAG_L001_R1_001

    :param filepath: a filename (possibly including a path)
    :type filepath: str
    :return: Unique sample ID
    :rtype: str
    """
    return os.path.basename(filepath).split('.fastq.gz')[0]


def get_sample_name_from_fastq_filename(filepath):
    """
    Takes:
    15-02380-CE11-T13-L1_AACCAG_L001_R1_001.fastq.gz

    Returns just the sample name:
    15-02380-CE11-T13-L1

    :param filepath: a filename (possibly including a path)
    :type filepath: str
    :return: Short sample name
    :rtype: str
    """
    parts = parse_sample_info_from_filename(filepath)
    return parts['sample_name']


def get_index_from_fastq_content(filepath):
    """
    Given a path to a FASTQ.gz file, returns the index sequence (assuming
    standard Illumina headers)

    :param filepath: The path to a gzipped FASTQ file.
    :type filepath: basestring
    :return: The index 'barcode' sequence.
    :rtype: str
    """
    with gzip.open(filepath, 'rb') as f:
        return f.readline().split(':').pop()
