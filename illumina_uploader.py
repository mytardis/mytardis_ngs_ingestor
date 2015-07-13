#!/usr/bin/env python
__author__ = 'Andrew Perry <Andrew.Perry@monash.edu.au>'

import sys
from urlparse import urlparse
import csv
from datetime import datetime

import os
from os.path import join, splitext
from dateutil import parser as dateparser
import subprocess

import xmltodict

from mytardis_uploader import MyTardisUploader
from mytardis_uploader import setup_logging, get_config, validate_config
# from mytardis_uploader import get_exclude_patterns_as_regex_list


def get_run_metadata(run_path):
    """

    :type run_path: str
    :rtype: dict
    """
    metadata = {}
    end_time, rta_version = rta_complete_parser(run_path)
    metadata['end_time'] = end_time.isoformat(' ')
    metadata['rta_version'] = rta_version
    metadata['run_path'] = run_path
    metadata['run_id'] = os.path.basename(run_path.rstrip(os.path.sep))

    runinfo_parameters = runinfo_parser(run_path)
    runinfo_parameters['rta_version'] = rta_version
    instrument_config = illumina_config_parser(run_path)

    metadata['instrument_model'] = instrument_config['instrument_model']
    metadata['title'] = "%s sequencing run %s" % (
        metadata['instrument_model'],
        metadata['run_id']
    )

    samplesheet, chemistry = get_samplesheet(join(run_path, 'SampleSheet.csv'))
    samplesheet_common_fields = {'chemistry': chemistry,
                                 'operator_name':
                                     samplesheet[0].get('Operator', '')
                                 }

    parameters = merge_dicts(runinfo_parameters,
                             instrument_config,
                             samplesheet_common_fields)
    parameter_list = dict_to_parameter_list(parameters)

    # parameter_list.append({u'name': param_name, u'value': param_value})
    schema = 'http://www.tardis.edu.au/schemas/ngs/illumina/run'
    metadata['parameter_sets'] = [{u'schema': schema,
                                   u'parameters': parameter_list}]

    return metadata


def get_project_metadata(proj_id,
                         run_metadata,
                         run_expt_link=None,
                         project_expt_link=None,
                         fastqc_summary_json='',
                         fastqc_version='',
                         schema='http://www.tardis.edu.au/schemas/ngs/project'):
    end_date = run_metadata['end_time'].split()[0]
    project_title = 'Sequencing Project, %s, %s' % (proj_id, end_date)
    if proj_id == 'Undetermined_indices':
        project_title = '%s, %s, %s' % (proj_id,
                                        run_metadata['run_id'],
                                        end_date)
    proj_metadata = {'title': project_title,
                     'description': project_title,
                     'institute': run_metadata['institute'],
                     'end_time': run_metadata['end_time'],
                     }

    # project_parameters = run_metadata
    # Project Dataset parameters are suffixed with '__project'
    # (to prevent name clashes with identical parameters at the
    #  Run Experiment level)
    # project_parameters = add_suffix_to_parameter_set(
    #    run_metadata['parameter_sets'][0]['parameters'], param_suffix)
    # project_parameter_list = dict_to_parameter_list(project_parameters)

    project_parameter_list = list(
        run_metadata['parameter_sets'][0]['parameters']
    )

    if run_expt_link is not None:
        project_parameter_list += dict_to_parameter_list(
            {'run_experiment': api_url_to_view_url(run_expt_link)}
        )
    if project_expt_link is not None:
        project_parameter_list += dict_to_parameter_list(
            {'project_experiment': api_url_to_view_url(project_expt_link)}
        )

    # This second (hidden) parameter_set, provides a summary of
    # FastQC results for every sample in the project, used for
    # rendering and overview table
    fastqc_summary_parameters = \
        dict_to_parameter_list({'fastqc_summary_json': fastqc_summary_json,
                                'fastqc_version': fastqc_version})

    fasqc_summary_parameterset = {
        u'schema':
        'http://www.tardis.edu.au/schemas/ngs/project/fastqc_summary',
        u'parameters': fastqc_summary_parameters
    }

    proj_metadata['parameter_sets'] = [{u'schema': schema,
                                        u'parameters': project_parameter_list},
                                       fasqc_summary_parameterset
                                       ]

    return proj_metadata

def api_url_to_view_url(apiurl):
    """
    Takes a MyTardis API URL of the form /v1/api/experiment/998 or
    and /v1/api/dataset/999, and returns the corresponding view URL,
    like /experiment/view/998 or /dataset/999.

    This is non-ideal and a bit of a hack, but there doesn't seem to
    be any easy way to achieve this with the existing MyTardis REST API
    (or without importing MyTardis server classes here to help).

    :type apiurl: str
    :rtype: str
    """
    pk = apiurl.strip('/').split('/')[-1:][0]
    if 'experiment' in apiurl:
        return '/experiment/view/%s' % pk
    if 'dataset' in apiurl:
        return '/dataset/%s' % pk

    raise ValueError('Not a valid MyTardis API URL: %s' % apiurl)

def add_suffix_to_parameter_set(parameters, suffix, divider='__'):
    """
    Adds a suffix ('__suffix') to the keys of a dictionary of MyTardis
    parameters. Returns a copy of the dict with suffixes on keys.
    (eg to prevent name clashes with identical parameters at the
        Run Experiment level)
    """
    suffixed_parameters = {}
    for p in parameters:
        suffixed_parameters[u'%s%s%s' % (p['name'],
                                         divider,
                                         suffix)] = p['value']
    return suffixed_parameters


def get_samplesheet(file_path):
    # FCID,Lane,SampleID,SampleRef,Index,Description,Control,Recipe,
    # Operator,SampleProject
    #
    # last line: #_IEMVERSION_3_TruSeq LT,,,,,,,,,
    with open(file_path, "rU") as f:
        reader = csv.DictReader(f)
        samples = [row for row in reader]
        chemistry = samples[-1:][0]['FCID'].split('_')[-1:][0]
        del samples[-1:]
    return samples, chemistry


def samplesheet_to_dict_by_id(samplesheet_rows):
    by_sample_id = {}
    for sample in samplesheet_rows:
        sample_id = sample['SampleID']
        by_sample_id[sample_id] = sample.copy()
        del by_sample_id[sample_id]['SampleID']

    return by_sample_id


def get_number_of_reads_fastq(filepath):
    """
    Count the number of reads in a (gzipped) FASTQ file.
    Assumes fours lines per read.

    :type filepath: str
    :rtype: int
    """

    # TODO: This should (optionally?) extract the number of reads from
    #       the associated FASTQC output - much faster !

    num = subprocess.check_output("zcat %s | echo $((`wc -l`/4))" % filepath,
                                  shell=True)
    return int(num.strip())


# Copypasta from: https://goo.gl/KpWo1w
# def unique(seq):
#     seen = set()
#     seen_add = seen.add
#     return [x for x in seq if not (x in seen or seen_add(x))]


def rta_complete_parser(run_path):
    """

    :type run_path: str
    :rtype: (datetime.DateTime, str)
    """
    # line in format:
    # 6/11/2014,20:00:49.935,Illumina RTA 1.17.20
    with open(join(run_path, "RTAComplete.txt"), 'r') as f:
        day, time, version = f.readline().split(',')

    date_time = dateparser.parse("%s %s" % (day, time))
    return date_time, version


def runinfo_parser(run_path):
    """

    Matches some or all of the fields defined in schema:
    http://www.tardis.edu.au/schemas/sequencing/illumina/run

    :type run_path: str
    :rtype: dict
    """
    with open(join(run_path, "RunInfo.xml"), 'r') as f:
        runinfo = xmltodict.parse(f)['RunInfo']['Run']

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

    Returns a dictionary with key/values of interest, where keys have
    the [section] from the Windows INI-style file are prepended,
    separated by colon. eg

    {"section_name:variable_name": "value"}

    :type run_path: str
    :rtype: dict
    """
    # we find the approriate config file
    config_filename = None
    for filename in os.listdir(join(run_path, 'Config')):
        if "_Effective.cfg" in filename:
            config_filename = filename
    if not config_filename:
        logger.error("Cannot find Config/*_Effective.cfg file")
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
    info = {}
    if 'system:instrumenttype' in allinfo:
        info['instrument_model'] = allinfo['system:instrumenttype']

    return info


def dict_to_parameter_list(d):
    """

    :type d: dict
    :rtype: list[dict]
    """
    return [{u'name': k, u'value': v} for k, v in d.items()]


# def merge_dicts(a, b):
#     merged = a.copy()
#     merged.update(b)
#     return merged

# http://stackoverflow.com/a/26853961
def merge_dicts(*dict_args):
    """
    Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.
    """
    result = {}
    for dictionary in dict_args:
        result.update(dictionary)
    return result


def create_run_experiment(metadata, uploader):
    """

    :type metadata: dict
    :type uploader: mytardis_uploader.MyTardisUploader
    :rtype: str
    """
    return uploader.create_experiment(
        metadata['title'],
        metadata['institute'],
        metadata['description'],
        end_time=metadata['end_time'],
        parameter_sets_list=metadata['parameter_sets'])


def create_project_experiment(metadata, uploader):
    """

    :type metadata: dict
    :type uploader: mytardis_uploader.MyTardisUploader
    :rtype: str
    """
    return uploader.create_experiment(metadata['title'],
                                      metadata['institute'],
                                      metadata['description'],
                                      end_time=metadata['end_time'],
                                      parameter_sets_list=metadata[
                                          'parameter_sets']
                                      )


def create_fastq_dataset(metadata, experiments, uploader):
    """

    :type metadata: dict
    :type experiments: list[str]
    :type uploader: mytardis_uploader.MyTardisUploader
    :rtype: str
    """
    return uploader.create_dataset(metadata['description'],
                                   experiments,
                                   parameter_sets_list=metadata[
                                       'parameter_sets']
                                   )


def register_project_fastq_datafiles(run_id,
                                     proj_path,
                                     samplesheet,
                                     dataset_url,
                                     uploader,
                                     fast_mode=False):

    schema = 'http://www.tardis.edu.au/schemas/ngs/fastq'

    sample_dict = samplesheet_to_dict_by_id(samplesheet)

    # Upload datafiles for the FASTQ reads in the project, for each Sample
    for sample_path, sample_id in get_sample_directories(proj_path):
        for fastq_path in get_fastq_read_files(sample_path):

                # sample_id may not be in the SampleSheet dict
                # if we are dealing the unaligned reads (eg sample_ids
                # lane1, lane2 etc)
                # So, we grab either the info for the sample_id, or
                # we use an empty dict (where all params will then be
                # the default missing values)

                # TODO: ensure we can also deal with unbarcoded runs,
                #       where Index is "NoIndex" and sample_ids are
                #       lane1, lane2 etc.

                sampleinfo = sample_dict.get(sample_id, {})

                reference_genome = sampleinfo.get('SampleRef', '')
                index_sequence = sampleinfo.get('Index', '')
                is_control = sampleinfo.get('Control', '')
                recipe = sampleinfo.get('Recipe', '')
                operator = sampleinfo.get('Operator', '')
                description = sampleinfo.get('Description', '')
                project = sampleinfo.get('SampleProject', '')

                if fast_mode:
                    number_of_reads = 0
                else:
                    number_of_reads = get_number_of_reads_fastq(fastq_path)

                parameters = {'run_id': run_id,
                              'sample_id': sample_id,
                              'reference_genome': reference_genome,
                              'index_sequence': index_sequence,
                              'is_control': is_control,
                              'recipe': recipe,
                              'operator': operator,
                              'description': description,
                              'project': project,
                              'number_of_reads': number_of_reads,
                              }
                datafile_parameter_sets = [
                    {u'schema': schema,
                     u'parameters': dict_to_parameter_list(parameters)}]

                replica_url = fastq_path
                if uploader.storage_mode == 'shared':
                    replica_url = get_shared_storage_replica_url(
                        uploader.storage_box_location,
                        fastq_path)

                if fast_mode:
                    md5_checksum = '__undetermined__'
                else:
                    md5_checksum = None  # will be calculated

                try:
                    uploader.upload_file(
                        fastq_path,
                        dataset_url,
                        parameter_sets_list=datafile_parameter_sets,
                        replica_url=replica_url,
                        md5_checksum=md5_checksum,
                    )
                except Exception, e:
                    logger.error("Failed to register Datafile: "
                                 "%s", fastq_path)
                    logger.debug("Exception: %s", e)
                    sys.exit(1)

                logger.info("Added Datafile: %s (%s)",
                            fastq_path,
                            dataset_url)


def get_sample_id_from_fastqc_filename(fastqc_zip_path):
    return os.path.basename(fastqc_zip_path).split('_')[0]


def register_project_fastqc_datafiles(run_id,
                                      proj_id,
                                      fastqc_out_dir,
                                      dataset_url,
                                      uploader,
                                      fast_mode=False):

    schema = 'http://www.tardis.edu.au/schemas/ngs/fastqc'

    # Upload datafiles for the FASTQC output files
    for fastqc_zip_path in get_fastqc_zip_files(fastqc_out_dir):

            sample_id = get_sample_id_from_fastqc_filename(fastqc_zip_path)
            fastqc_data_tables = parse_fastqc_data_txt(fastqc_zip_path)
            parameters = {'run_id': run_id,
                          'project': proj_id,
                          'sample_id': sample_id,
                          'fastqc_version': fastqc_data_tables['version']}

            datafile_parameter_sets = [
                {u'schema': schema,
                 u'parameters': dict_to_parameter_list(parameters)}]

            replica_url = fastqc_zip_path
            if uploader.storage_mode == 'shared':
                replica_url = get_shared_storage_replica_url(
                    uploader.storage_box_location,
                    fastqc_zip_path)

            if fast_mode:
                md5_checksum = '__undetermined__'
            else:
                md5_checksum = None  # will be calculated

            try:
                uploader.upload_file(
                    fastqc_zip_path,
                    dataset_url,
                    parameter_sets_list=datafile_parameter_sets,
                    replica_url=replica_url,
                    md5_checksum=md5_checksum,
                )
            except Exception, e:
                logger.error("Failed to register Datafile: "
                             "%s", fastqc_zip_path)
                logger.debug("Exception: %s", e)
                sys.exit(1)

            logger.info("Added Datafile: %s (%s)",
                        fastqc_zip_path,
                        dataset_url)


def upload_fastqc_reports(fastqc_out_dir, dataset_url, options):

    for fastqc_zip_path in get_fastqc_zip_files(fastqc_out_dir):

        # Since the zipped FastQC output may be migrated to a storage location
        # on tape (with delayed access), we need to host 'live' html reports
        # at a location that will always be immediately available

        live_box_uploader = MyTardisUploader(
            options.url,
            options.username,
            options.password,
            storage_mode='upload',
            # storage_box_location='/data/cached',
            storage_box_name='live')

        fastqc_data_tables = parse_fastqc_data_txt(fastqc_zip_path)
        fqc_version = fastqc_data_tables['version']

        # Depending on the version of FastQC used and specific
        #       commandline options, we may or may not have an HTML
        #       report with inline (Base64) images. As such:
        #       * if html file exists & FastQC version is,
        #         >= 0.11.3, upload that file (assume inline images)
        #       * if FastQC version is < 0.11.3, extra the zip file to /tmp,
        #         create an inline images version, then upload that
        #
        #       Currently FastQC always generates a zip file alongside
        #       any other output, irrespective of command line options.
        #       For this reason, we always just extract from the zip, since
        #       we know it should be there
        with TmpZipExtract(fastqc_zip_path) as tmp_path:
            fqc_id = splitext(os.path.basename(fastqc_zip_path))[0]
            report_dir = join(tmp_path, fqc_id)
            report_file = join(report_dir, 'fastqc_report.html')
            from distutils.version import LooseVersion as SoftwareVersion
            inline_report_filename = '%s.html' % fqc_id
            inline_report_abspath = join(report_dir, inline_report_filename)
            if SoftwareVersion(fqc_version) < SoftwareVersion('0.11.3'):
                # convert fastqc_report.html to version with inline images
                from standalone_html import make_html_images_inline
                make_html_images_inline(report_file, inline_report_abspath)
            else:
                os.rename(report_file, inline_report_abspath)

            location = live_box_uploader.upload_file(inline_report_abspath,
                                                     dataset_url)


def get_fastqc_summary_table_for_project(fastqc_out_dir):

    project_summary = []
    # Upload datafiles for the FASTQC output files
    for fastqc_zip_path in get_fastqc_zip_files(fastqc_out_dir):
        sample_id = get_sample_id_from_fastqc_filename(fastqc_zip_path)
        qc_pass_fail_table = parse_fastqc_summary_txt(fastqc_zip_path)

        # TODO: Currently unused, except for extracting the FastQC version
        qc_data_table = parse_fastqc_data_txt(fastqc_zip_path)

        project_summary.append((sample_id, qc_pass_fail_table))

    return project_summary, qc_data_table['version']

# def get_project_directories(bcl2fastq_out_dir):
#     """
#     Returns an iterator that gives "Project_" directories from
#     the bcl2fastq output directory.
#
#     :type bcl2fastq_out_dir: str
#     :return: Iterator
#     """
#     for item in os.listdir(bcl2fastq_out_dir):
#         proj_path = os.path.join(bcl2fastq_out_dir, item)
#         # we look for directories named Project_*
#         # NOTE: an alternative method would be to use the projects explicitly
#         # listed in SameleSheet.csv, treating it as the one source of
#         # truth and avoiding cases where someone decided to create
#         # an additional Project_ directory (eg Project_Bla.old)
#         if os.path.isdir(proj_path) and ('Project_' in item):
#             yield proj_path


def get_bcl2fastq_output_dir(options, run_metadata):
    return options.bcl2fastq_output_path.format(
        run_id=run_metadata['run_id'],
        run_path=run_metadata['run_path']
    )


def get_sample_directories(project_path):
    """
    Returns an iterator that gives tuples of ("Sample_" directory, Sample ID)
    from within a bcl2fastq Project_ directory.

    :type project_path: str
    :return: Iterator
    """
    for item in os.listdir(project_path):
        sample_path = os.path.join(project_path, item)
        if os.path.isdir(sample_path) and ('Sample_' in item):
            yield sample_path, item.split('Sample_')[1]


def get_fastq_read_files(sample_path):
    """
    Returns an iterator that gives gzipped FASTQ read files from a
    a bcl2fastq Project_*/Sample_* directory.

    :type sample_path: str
    :return: Iterator
    """
    for item in os.listdir(sample_path):
        fastq_path = os.path.join(sample_path, item)

        # NOTE: an alternative would be to use SampleSheet.csv
        # to construct the expected output filenames
        # Illumina FASTQ files use the following naming scheme:
        # <sample id>_<index>_L<lane>_R<read number>_<setnumber>.fastq.gz
        # For example, the following is a valid FASTQ file name:
        # NA10831_ATCACG_L002_R1_001.fastq.gz
        # Note that lane and set numbers are 0-padded to 3 digits.
        # In the case of non-multiplexed runs, <sample name> will be replaced
        # with the lane numbers (lane1, lane2, ..., lane8) and <index> will be
        #  replaced with "NoIndex".
        # read_number = 1  # or 2?
        # set_number = 1
        # fn = '%s_%s_L%03d_R%d_%03d' % (ss['SampleID'],
        #                                ss['Index'],
        #                                ss['Lane'],
        #                                read_number,
        #                                set_number)

        # we just return things with the .gz extension
        if os.path.splitext(item)[-1:][0] == '.gz' and \
           os.path.isfile(fastq_path):
            yield fastq_path


def get_fastqc_zip_files(fastqc_out_path):
    """
    Returns an iterator that gives zipped FASTQC result files.

    :type fastqc_out_path: str
    :return: Iterator
    """
    for item in os.listdir(fastqc_out_path):
        fastqc_path = os.path.join(fastqc_out_path, item)

        # we just return things with the .zip extension
        if os.path.splitext(item)[-1:][0] == '.zip' and \
           os.path.isfile(fastqc_path):
            yield fastqc_path


def get_fastqc_output_directory(proj_path):
    return join(proj_path, 'FastQC.out')


def run_fastqc(fastq_paths,
               output_directory=None,
               fastqc_bin=None,
               threads=2,
               extra_options=''):

    tmp_dir = None
    if not output_directory:
        try:
            from tempfile import mkdtemp
            tmp_dir = mkdtemp()
            output_directory = tmp_dir
        except OSError:
            logger.error('FastQC - failed to create temp directory.')
            return None

    if not fastqc_bin:
        # We will assume it's on path with the default name
        fastqc_bin = 'fastqc'

    cmd = '%s %s --noextract --threads %d --outdir %s %s' % \
          (fastqc_bin,
           extra_options,
           threads,
           output_directory,
           ' '.join(fastq_paths))

    logger.info('Running FastQC on: %s', ', '.join(fastq_paths))

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
            import shutil
            shutil.rmtree(output_directory)
            return None

    return output_directory


def run_fastqc_on_project(proj_path,
                          output_directory=None,
                          fastqc_bin=None,
                          threads=2):

    if output_directory is None:
        output_directory = get_fastqc_output_directory(proj_path)
        if os.path.exists(output_directory):
            logger.error("FastQC - output directory already exists: %s",
                         output_directory)
            return None
        try:
            os.mkdir(output_directory)
        except OSError:
            logger.error("FastQC - couldn't create output directory: %s",
                         output_directory)
            return None

    if (not os.path.exists(output_directory)) or \
       (not os.path.isdir(output_directory)):
        logger.error("FastQC - output path %s isn't a directory "
                     "or doesn't exist.",
                     output_directory)
        return None

    fastq_files = []
    for sample_path, sample_id in get_sample_directories(proj_path):
        fastq_files.extend([f for f in get_fastq_read_files(sample_path)])

    fastqc_out = run_fastqc(fastq_files,
                            output_directory=output_directory,
                            fastqc_bin=fastqc_bin,
                            threads=threads)
    return fastqc_out


def file_from_zip(zip_file_path, filename, mode='r'):
    # eg rootpath =
    # FastQC.out/15-02380-CE11-T13-L1_AACCAG_L001_R1_001_fastqc.zip
    # ie   fastq_filename + '_fastqc.zip'

    from fs.opener import opener
    # from fs.zipfs import ZipFS

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


def parse_fastqc_summary_txt(zip_file_path):
    """
    Extract the overall PASS/WARN/FAIL summary.txt table from FastQC results
    contained in a zip file.
    """
    def parse(fh):
        summary = [tuple(line.strip().split('\t')) for line in fh]
        return summary

    return parse_file_from_zip(zip_file_path,
                               'summary.txt',
                               parse)


def parse_fastqc_data_txt(zip_file_path):
    """
    Extract the tables in fastqc_data.txt from FastQC results contained
    in a zip file.
    """
    def parse(fh):
        data = {'version': fh.readline().split('\t')[1].strip()}
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

    return parse_file_from_zip(zip_file_path,
                               'fastqc_data.txt',
                               parse)


class TmpZipExtract:
    """
    A Context Manager class that unzips a zip file to a temporary path
    and cleans up afterwards, when used with the 'with' statement.
    """
    def __init__(self, zip_file_path):
        self.zip_file_path = zip_file_path
        self.tmpdir = None

    def __enter__(self):
        from tempfile import mkdtemp
        tmpdir = mkdtemp()
        from zipfile import ZipFile
        zf = ZipFile(self.zip_file_path)
        zf.extractall(tmpdir)
        self.tmpdir = tmpdir
        return tmpdir

    def __exit__(self, type, value, traceback):
        import shutil
        shutil.rmtree(self.tmpdir)

def extract_fastqc_zip(zip_file_path):
    """
    Extract a zip file to a temporary path. Returns the path.

    :type zip_file_path: str
    :rtype: str
    """
    from tempfile import mkdtemp
    tmpdir = mkdtemp()
    from zipfile import ZipFile
    zf = ZipFile(zip_file_path)
    zf.extractall(tmpdir)
    return tmpdir

    #file_from_zip(zip_file_path, 'fastqc_report.html')


def get_shared_storage_replica_url(storage_box_location, file_path):
    """
    Generates a 'replica_url' for a storage box location when
    using 'shared' storage mode.

    eg, if storage box base path is:
        /data/bigstorage/
    and absolute file path is
        /data/bigstorage/expt1/dataset1/file.txt
    then replica_url will be:
        expt1/dataset1/file.txt

    :type file_path: str
    :rtype: str
    """
    if storage_box_location:
        replica_url = os.path.relpath(file_path,
                                      storage_box_location)
        return replica_url


def run_main():
    logger = setup_logging()
    global logger

    def extra_config_options(argparser):
        argparser.add_argument('--threads',
                               dest='threads',
                               type=int,
                               metavar='THREADS')
        argparser.add_argument('--run-fastqc',
                               dest='run_fastqc',
                               type=bool,
                               default=False,
                               metavar='RUN_FASTQC')
        argparser.add_argument('--fastqc-bin',
                               dest='fastqc_bin',
                               type=str,
                               metavar='FASTQC_BIN')
        argparser.add_argument('--bcl2fastq-output-path',
                               dest='bcl2fastq_output_path',
                               default='{run_path}/{run_id}.bcl2fastq',
                               type=str,
                               metavar='BCL2FASTQ_OUTPUT_PATH')

    parser, options = get_config(add_extra_options_fn=extra_config_options)

    validate_config(parser, options)

    # exclude_patterns = \
    #     get_exclude_patterns_as_regex_list(options.exclude)

    uploader = MyTardisUploader(
        options.url,
        options.username,
        options.password,
        storage_mode=options.storage_mode,
        storage_box_location=options.storage_base_path
    )

    # Create an Experiment representing the overall sequencing run
    run_path = options.path
    run_metadata = get_run_metadata(run_path)

    run_metadata['institute'] = options.institute
    run_metadata['description'] = options.description
    if not run_metadata['description']:
        run_metadata['description'] = "Automatically ingested by %s on %s" % \
                                      (uploader.user_agent,
                                       datetime.now().isoformat(' '))

    try:
        run_expt_url = create_run_experiment(run_metadata, uploader)
    except Exception, e:
        logger.error("Failed to create Experiment for sequencing run: %s",
                     run_path)
        logger.error("Exception: %s: %s", type(e).__name__, e)
        sys.exit(1)

    # Take just the path of the experiment, eg /api/v1/experiment/187/
    run_expt_url = urlparse(run_expt_url).path

    logger.info("Created run experiment: %s (%s)",
                run_metadata['run_id'],
                run_expt_url)

    samplesheet, chemistry = get_samplesheet(join(run_path, 'SampleSheet.csv'))

    # Get unique project names
    projects = list(set([s['SampleProject'] for s in samplesheet]))
    # Treat the 'Undetermined_indices' directory as a (special) project
    projects.append('Undetermined_indices')

    for proj_id in projects:
        proj_metadata = get_project_metadata(
            proj_id,
            run_metadata,
            run_expt_link=run_expt_url,
            schema='http://www.tardis.edu.au/schemas/ngs/project')

        # The directory where bcl2fastq puts its output,
        # in Project_* directories
        bcl2fastq_output_dir = get_bcl2fastq_output_dir(options,
                                                        run_metadata)

        proj_path = os.path.join(bcl2fastq_output_dir, 'Project_' + proj_id)
        fastqc_out_dir = get_fastqc_output_directory(proj_path)

        if proj_id == 'Undetermined_indices':
            proj_path = os.path.join(bcl2fastq_output_dir, proj_id)

        # run FastQC if output doesn't exist
        if options.run_fastqc and not options.fast:
            run_fastqc_on_project(proj_path,
                                  fastqc_bin=options.fastqc_bin,
                                  threads=int(options.threads))

        fqc_summary = []
        if proj_id != 'Undetermined_indices':
            fqc_summary, fqc_version = \
                get_fastqc_summary_table_for_project(fastqc_out_dir)

        import json
        fqc_summary_json = json.dumps(fqc_summary)

        # Create an Experiment for each real Project in the run
        # (except Undetermined_indices Datasets which will only be associated
        #  with the parent 'run' Experiment)
        if proj_id != 'Undetermined_indices':
            try:
                project_url = create_project_experiment(proj_metadata, uploader)
            except Exception, e:
                logger.error("Failed to create Experiment for project: %s",
                             proj_id)
                logger.debug("Exception: %s", e)
                sys.exit(1)

            project_url = urlparse(project_url).path

            logger.info("Created Project Experiment: %s (%s)",
                        project_url,
                        proj_id)

        # samples_in_project = [s['SampleID'] for s in samplesheet]
        # sample_desc = 'FASTQ reads for samples '+', '.join(samples_in_project)

        dataset_metadata = get_project_metadata(
            proj_id,
            run_metadata,
            run_expt_link=run_expt_url,
            project_expt_link=project_url,
            fastqc_summary_json=fqc_summary_json,
            fastqc_version=fqc_version,
            schema='http://www.tardis.edu.au/schemas/ngs/raw_reads')

        # Create a Dataset for each project, associated with both
        # the overall run Experiment, and the project Experiment
        try:
            # We associate the Undetermined_indices FASTQ reads only
            # with the overall 'run' Experiment (unlike proper Projects
            # which also have their own Project Experiment).
            if proj_id == 'Undetermined_indices':
                parent_expt_urls = [run_expt_url]
            else:
                parent_expt_urls = [project_url, run_expt_url]

            dataset_url = create_fastq_dataset(dataset_metadata,
                                               parent_expt_urls,
                                               uploader)
        except Exception, e:
            logger.error("Failed to create Dataset for project: %s",
                         proj_id)
            logger.debug("Exception: %s", e)
            sys.exit(1)

        # Take just the path, eg: /api/v1/dataset/363
        dataset_url = urlparse(dataset_url).path

        logger.info("Created FASTQ dataset: %s (%s)", dataset_url, proj_id)

        register_project_fastq_datafiles(run_metadata['run_id'],
                                         proj_path,
                                         samplesheet,
                                         dataset_url,
                                         uploader,
                                         fast_mode=options.fast)

        if proj_id != 'Undetermined_indices':
            register_project_fastqc_datafiles(run_metadata['run_id'],
                                              proj_id,
                                              fastqc_out_dir,
                                              dataset_url,
                                              uploader,
                                              fast_mode=options.fast)

            upload_fastqc_reports(fastqc_out_dir, dataset_url, options)

    logger.info("Ingestion of run %s complete !", run_metadata['run_id'])
    sys.exit(0)

if __name__ == "__main__":
    run_main()
