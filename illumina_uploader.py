#!/usr/bin/env python
__author__ = 'Andrew Perry <Andrew.Perry@monash.edu.au>'

import sys
from urlparse import urlparse
import csv
from datetime import datetime

import os
from os.path import join
from dateutil import parser as dateparser

import xmltodict

from mytardis_uploader import MyTardisUploader
from mytardis_uploader import setup_logging, get_config, validate_config
from mytardis_uploader import get_exclude_patterns_as_regex_list


def get_run_metadata(run_path):
    """

    :param run_path: str
    :rtype: dict
    """
    metadata = {}
    end_time, rta_version = rta_complete_parser(run_path)
    metadata['end_time'] = end_time.isoformat(' ')
    metadata['rta_version'] = rta_version
    metadata['run_dir'] = os.path.relpath(
        options.path,
        options.storage_base_path
    )

    runinfo_parameters = runinfo_parser(run_path)
    instrument_config = illumina_config_parser(run_path)

    metadata['instrument_model'] = instrument_config['instrument_model']
    metadata['title'] = "%s sequencing run %s" % (
        metadata['instrument_model'],
        metadata['run_dir'])

    parameters = merge_dicts(runinfo_parameters, instrument_config)
    parameters = dict_to_parameter_list(parameters)

    # parameter_list.append({u'name': param_name, u'value': param_value})
    schema = 'http://www.tardis.edu.au/schemas/sequencing/illumina/run'
    metadata['parameter_sets'] = [{u'schema': schema,
                                  u'parameters': parameters}]

    return metadata

def get_samplesheet(file_path):
    # FCID,Lane,SampleID,SampleRef,Index,Description,Control,Recipe,
    # Operator,SampleProject
    #
    # last line: #_IEMVERSION_3_TruSeq LT,,,,,,,,,
    with open(file_path, "rU") as f:
        reader = csv.DictReader(f)
        samples = [row for row in reader]
        chemistry = samples[-1:][0]['FCID'].split('_')[-1:]
        del samples[-1:]
    return samples, chemistry

# Copypasta from: https://goo.gl/KpWo1w
# def unique(seq):
#     seen = set()
#     seen_add = seen.add
#     return [x for x in seq if not (x in seen or seen_add(x))]

def get_samplesheet_projects(samplesheet):
    """

    :param samplesheet: dict
    :rtype: list[str]
    """
    for sample in samplesheet:
        pass


def rta_complete_parser(run_path):
    """

    :param run_path: str
    :rtype: (datetime.DateTime, str)
    """
    # line in format:
    # 6/11/2014,20:00:49.935,Illumina RTA 1.17.20
    with open(join(run_path, "RTAComplete.txt"), 'r') as f:
        day, time, version = f.readline().split(',')

    datetime = dateparser.parse("%s %s" % (day, time))
    return datetime, version

def runinfo_parser(run_path):
    """

    Matches some or all of the fields defined in schema:
    http://www.tardis.edu.au/schemas/sequencing/illumina/run

    :param run_path:
    :rtype: dict
    """
    with open(join(run_path, "RunInfo.xml"), 'r') as f:
        runinfo = xmltodict.parse(f)['RunInfo']['Run']

    info = {
            u'run_id': runinfo['@Id'],
            u'run_number': runinfo['@Number'],
            u'flowcell_id': runinfo['Flowcell'],
            u'instrument_id': runinfo['Instrument']
           }

    reads = runinfo['Reads']['Read']

    # NOTE: We don't do this since generating many
    #       numbered parameter names isn't cleanly
    #       compatible with the rigidly defined parameter
    #       sets defined by a MyTardis schema.
    # for read in reads:
    #     readtype = u'read_'
    #     if read['@IsIndexedRead'] == 'Y':
    #         readtype = u'read_index_'
    #     readtype += read['@Number'] + u'_cycles'
    #     info.append({readtype: read['@NumCycles']})

    cycle_list = []
    index_reads = []
    for read in reads:
        if read['@IsIndexedRead'] == 'Y':
            index_reads.append(read['@Number'])
        # we assume reads are always in order of Number
        cycle_list.append(read['@NumCycles'])

    info['read_cycles'] = ', '.join(cycle_list)
    info['index_reads'] = ', '.join(index_reads)

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

    :param run_path: str
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

    :param d: list[dict]
    :rtype: list[dict]
    """
    return [{u'name': k, u'value': v} for k, v in d.items()]

def merge_dicts(a, b):
    merged = a.copy()
    merged.update(b)
    return merged

def create_run_experiment(metadata, uploader):
    """

    :param metadata: dict
    :param uploader: mytardis_uploader.MyTardisUploader
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

    :param metadata: dict
    :param uploader: mytardis_uploader.MyTardisUploader
    :rtype: str
    """
    return uploader.create_experiment(metadata['title'],
                                      metadata['institute'],
                                      metadata['description'],
                                      end_time=metadata['end_time'])

def create_fastq_dataset(metadata, experiments, uploader):
    """

    :param metadata: dict
    :param experiments: list[str]
    :param uploader: mytardis_uploader.MyTardisUploader
    :rtype: str
    """
    return uploader.create_dataset(metadata['description'],
                                   experiments,
                                   parameter_sets_list=None)

def get_project_directories(bcl2fastq_out_dir):
    """
    Returns an iterator that gives "Project_" directories from
    the bcl2fastq output directory.

    :param bcl2fastq_out_dir: str
    :return: Iterator
    """
    for item in os.listdir(bcl2fastq_out_dir):
        proj_path = os.path.join(bcl2fastq_out_dir, item)
        # we look for directories named Project_*
        # NOTE: an alternative method would be to use the projects explicitly
        # listed in SameleSheet.csv, treating it as the one source of
        # truth and avoiding cases where someone decided to create
        # an additional Project_ directory (eg Project_Bla.old)
        if os.path.isdir(proj_path) and ('Project_' in item):
            yield proj_path

def get_sample_directories(project_path):
    """
    Returns an iterator that gives "Sample_" directories from within a
    a bcl2fastq Project_ directory.

    :param project_path: str
    :return: Iterator
    """
    for item in os.listdir(project_path):
        sample_path = os.path.join(proj_path, item)
        if os.path.isdir(sample_path) and ('Sample_' in item):
            yield sample_path

def get_fastq_read_files(sample_path):
    """
    Returns an iterator that gives gzipped FASTQ read files from a
    a bcl2fastq Project_*/Sample_* directory.

    :param sample_path: str
    :return: Iterator
    """
    for item in os.listdir(sample_path):
        fastq_path = os.path.join(sample_path, item)
        # we just return things with the .gz extension
        # NOTE: an alternative would be to use SampleSheet.csv
        # to construct the expected output filenames
        if os.path.splitext(item)[-1:][0] == '.gz' and \
           os.path.isfile(fastq_path):
            yield fastq_path

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

    :param uploader: str
    :param file_path: str
    :rtype: str
    """
    if storage_box_location:
        replica_url = os.path.relpath(file_path,
                                      storage_box_location)
        return replica_url

if __name__ == "__main__":
    logger = setup_logging()

    parser, options = get_config()

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

    # create an Experiment representing the overall sequencing run
    run_path = options.path
    run_metadata = get_run_metadata(run_path)
    run_metadata['institute'] = options.institute
    run_metadata['description'] = options.description
    if not run_metadata['description']:
        run_metadata['description'] = "Automatically ingested by %s on %s" % \
                           (uploader.user_agent, datetime.now().isoformat(' '))

    try:
        run_expt_url = create_run_experiment(run_metadata, uploader)
    except Exception, e:
        logger.error("Failed to create Experiment for sequencing run: %s",
                     run_path)
        logger.error("Exception: %s", e)
        sys.exit(1)

    # take just the path of the experiment, eg /api/v1/experiment/187/
    run_expt_url = urlparse(run_expt_url).path

    logger.info("Created run experiment: %s (%s)",
                run_metadata['run_dir'],
                run_expt_url)

    samplesheet, chemistry = get_samplesheet(join(run_path, 'SampleSheet.csv'))

    # get unique project names
    projects = list(set([s['SampleProject'] for s in samplesheet]))

    for project in projects:
        # create an Experiment for each project in the run
        end_date = run_metadata['end_time'].split()[0]
        project_title = 'Sequencing Project, %s, %s' % (project, end_date)
        metadata = {'title': project_title,
                    'description': project_title,
                    'institute': run_metadata['institute'],
                    'end_time': run_metadata['end_time'],
                    }

        try:
            project_url = create_project_experiment(metadata, uploader)
        except Exception, e:
            logger.error("Failed to create Experiment for project: %s",
                         project)
            logger.debug("Exception: %s", e)
            sys.exit(1)

        project_url = urlparse(project_url).path

        logger.info("Created Project Experiment: %s (%s)", project_url, project)

        # create a Dataset for each project, associated with both
        # the overall run Experiment, and the project Experiment
        samples_in_project = [s['SampleID'] for s in samplesheet]
        sample_desc = 'FASTQ reads for samples ' + ', '.join(samples_in_project)
        metadata = {'description': project_title}

        try:
            dataset_url = create_fastq_dataset(metadata,
                                               [project_url, run_expt_url],
                                               uploader)
        except Exception, e:
            logger.error("Failed to create Dataset for project: %s",
                         project)
            logger.debug("Exception: %s", e)
            sys.exit(1)

        dataset_url = urlparse(dataset_url).path

        logger.info("Created FASTQ dataset: %s (%s)", dataset_url, project)

        # the directory where bcl2fastq puts its output,
        # in Project_* directories
        bcl2fastq_output_dir = join(run_path, run_metadata['run_dir']+'.pc')

        # upload datafiles for the FASTQ reads in the project,
        # for each Sample
        for proj_path in get_project_directories(bcl2fastq_output_dir):
            for sample_path in get_sample_directories(proj_path):
                for fastq_path in get_fastq_read_files(sample_path):
                        replica_url = fastq_path
                        if uploader.storage_mode == 'shared':
                            replica_url = get_shared_storage_replica_url(
                                uploader.storage_box_location,
                                fastq_path)

                        try:
                            uploader.upload_file(
                                fastq_path, dataset_url,
                                parameter_sets_list=None,
                                replica_url=replica_url)
                        except Exception, e:
                            logger.error("Failed to register Datafile: "
                                         "%s", fastq_path)
                            logger.debug("Exception: %s", e)
                            sys.exit(1)

                        logger.info("Added Datafile: %s (%s)",
                                    fastq_path,
                                    dataset_url)

    logger.info("Ingestion of run %s complete !", run_metadata['run_dir'])
    sys.exit(0)