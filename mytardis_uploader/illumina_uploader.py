#!/usr/bin/env python
__author__ = 'Andrew Perry <Andrew.Perry@monash.edu.au>'

import sys
import shutil
from urlparse import urlparse
import csv
from datetime import datetime
import subprocess
from tempfile import mkdtemp
import atexit
from urlparse import urljoin

import os
from os.path import join, splitext, exists, isdir, isfile
from dateutil import parser as dateparser
import xmltodict
import json

from semantic_version import Version as SemanticVersion
from distutils.version import LooseVersion

import mytardis_uploader
from mytardis_uploader import MyTardisUploader
from mytardis_uploader import setup_logging, get_config, validate_config
# from mytardis_uploader import get_exclude_patterns_as_regex_list

import models
from mytardis_models import Experiment, Dataset, DataFile
from models import DemultiplexedSamplesBase, FastqcOutputBase, \
    FastqcReportsBase, FastqRawReadsBase, HiddenFastqcProjectSummaryBase, \
    IlluminaSequencingRunBase, NucleotideRawReadsDatasetBase, \
    IlluminaRunConfigBase

# a module level list of temporary directories that have been
# created, so these can be cleaned up upon premature exit
TMPDIRS = []


class DemultiplexedSamples(DemultiplexedSamplesBase):
    pass


class FastqcOutput(FastqcOutputBase):
    pass


class FastqcReports(FastqcReportsBase):
    pass


class FastqRawReads(FastqRawReadsBase):
    pass


class HiddenFastqcProjectSummary(HiddenFastqcProjectSummaryBase):
    pass


class IlluminaSequencingRun(IlluminaSequencingRunBase):
    pass


class NucleotideRawReadsDataset(NucleotideRawReadsDatasetBase):
    pass


class IlluminaRunConfig(IlluminaRunConfigBase):
    pass


def _format_day(dt):
    return dt.strftime("%d-%b-%Y")


def create_run_experiment_object(run_path):
    """

    :type run_path: str
    :rtype: Experiment
    """

    # TODO: grab start time metadata from somewhere
    #       (eg maybe from Logs/CycleTimes.txt)
    end_time, rta_version = rta_complete_parser(run_path)

    runinfo_parameters = runinfo_parser(run_path)
    instrument_config = illumina_config_parser(run_path)
    instrument_model = instrument_config.get('system:instrumenttype', '')
    instrument_id = runinfo_parameters.get('instrument_id', '')
    samplesheet, chemistry = parse_samplesheet(join(run_path,
                                                    'SampleSheet.csv'))

    # the MyTardis ParameterSet
    run = IlluminaSequencingRun()
    run.from_dict(runinfo_parameters)
    run.instrument_model = instrument_model
    run.instrument_id = instrument_id
    # run.run_path = run_path
    run.run_id = get_run_id_from_path(run_path)
    run.rta_version = rta_version
    run.chemistry = chemistry
    run.operator_name = samplesheet[0].get('Operator', '')
    run._samplesheet = samplesheet

    # the MyTardis Experiment
    expt = Experiment()
    expt.title = "%s sequencing run %s" % (run.instrument_model, run.run_id)
    expt.end_time = end_time
    expt.parameters = run
    expt.run_path = run_path

    return expt


def create_project_experiment_object(
        proj_id,
        run_expt,
        run_expt_link=None):

    """

    :type proj_id: str
    :type run_expt: Experiment
    :type run_expt_link: str
    :rtype: Experiment
    """

    end_date = _format_day(run_expt.end_time)
    project_title = 'Sequencing Project, %s, %s' % (proj_id, end_date)
    if proj_id == 'Undetermined_indices':
        project_title = '%s, %s, %s' % (proj_id,
                                        run_expt.parameters.run_id,
                                        end_date)

    proj_expt = Experiment()
    proj_expt.title = project_title
    proj_expt.description = project_title
    proj_expt.institution_name = run_expt.institution_name
    proj_expt.end_time = run_expt.end_time
    proj_expt.run_path = run_expt.run_path

    project_params = DemultiplexedSamples()
    project_params.from_dict(run_expt.parameters.to_dict(),
                             existing_only=False)
    project_params.project_id = proj_id
    if run_expt_link is not None:
        project_params.run_experiment = run_expt_link
    proj_expt.parameter_sets.append(project_params)

    return proj_expt


def create_fastqc_dataset_object(run_id,
                                 proj_id,
                                 end_time,
                                 experiments,
                                 raw_reads_dataset_url,
                                 fastqc_summary=None):
    end_date = _format_day(end_time)
    fqc_dataset_title = "FastQC reports for Project %s, %s" % \
                        (proj_id, end_date)
    if proj_id == 'Undetermined_indices':
        fqc_dataset_title = "FastQC reports for %s, %s" % (proj_id, end_date)

    fqc_dataset = Dataset()
    fqc_dataset.experiments = experiments
    fqc_dataset.description = fqc_dataset_title

    fqc_version = ''
    if fastqc_summary is not None:
        fqc_version = fastqc_summary.get('fastqc_version', '')

    params = {
        'run_id': run_id,
        'project': proj_id,
        'raw_reads_dataset': raw_reads_dataset_url,
        'fastqc_version': fqc_version,
    }
    fqc_params = FastqcReports()
    fqc_params.from_dict(params)

    fqc_dataset.parameter_sets.append(fqc_params)

    # This second (hidden) parameter_set, provides a summary of
    # FastQC results for every sample in the project, used for
    # rendering and overview table
    if fastqc_summary is not None:
        fastqc_summary_params = HiddenFastqcProjectSummary()
        fastqc_summary_params.hidden_fastqc_summary_json = \
            json.dumps(fastqc_summary)
        fastqc_summary_params.fastqc_version = fqc_version
        fqc_dataset.parameter_sets.append(fastqc_summary_params)

    return fqc_dataset


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


def parse_samplesheet(file_path):
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


def filter_samplesheet_by_project(file_path, proj_id,
                                  project_column_label='SampleProject'):
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
        header = f.readline().strip()
        s = header.split(',')
        project_column_index = s.index(project_column_label)
        outlines.append(header+'\r\n')
        for l in f:
            s = l.strip().split(',')
            if s[project_column_index] == proj_id or l[0] == '#':
                outlines.append(l)
    return outlines


def samplesheet_to_dict_by_id(samplesheet_rows):
    by_sample_id = {}
    for sample in samplesheet_rows:
        sample_id = sample['SampleID']
        by_sample_id[sample_id] = sample.copy()
        del by_sample_id[sample_id]['SampleID']

    return by_sample_id


def get_project_ids_from_samplesheet(samplesheet,
                                     include_no_index_name=None):

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
    num = subprocess.check_output("zcat -f %s | "
                                  "head -n 2  | "
                                  "tail -n 1  | "
                                  "wc --chars" % filepath,
                                  shell=True)
    return int(num.strip()) - 1

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

    # select just the keys of interest to return
    # info = {}
    # if 'system:instrumenttype' in allinfo:
    #    info['instrument_model'] = allinfo['system:instrumenttype']

    return allinfo


def dict_to_parameter_list(d):
    """

    :type d: dict
    :rtype: list[dict]
    """
    return [{u'name': k, u'value': v} for k, v in d.items()]


def dict_to_parameter_set(d, schema):
    return {u'schema': schema,
            u'parameters': dict_to_parameter_list(d)}


def parameter_set_to_dict(param_set):
    d = {}
    for param in param_set['parameters']:
        d[param['name']] = param['value']

    return d


def query_experiment_parameterset(self,
                                  uploader,
                                  schema_namespace=None,
                                  parameter_name=None,
                                  value=None):
    query_params = {}
    if isinstance(schema_namespace, basestring):
        query_params[u'name__schema__namespace'] = schema_namespace
    if isinstance(parameter_name, basestring):
        query_params[u'name__name'] = parameter_name
    if isinstance(value, (float, int)):
        query_params[u'numeric_value'] = value
    if isinstance(value, (basestring)):
        query_params[u'string_value'] = value
    if isinstance(value, datetime.datetime):
        query_params[u'datetime_value'] = value.isoformat()

    response = uploader.do_get_request('experimentparameter', query_params)
    return response.json()


def get_experiments_from_server_by_run_id(uploader, run_id,
                                          schema_namespace):
    """
    Query the server for all Experiments (run and project) with a matching
    run_id Parameter. Returns a list of integer object IDs.

    :param uploader: An uploader object instance
    :type uploader: MyTardisUploader
    :param run_id: The unique run ID (eg 150225_SNL177_0111_AHFNLKADXX)
    :type run_id: basestring
    :param schema_namespace: The namespace of the MyTardis schema (a URL)
    :type schema_namespace: basestring
    :return: Object IDs of all matching experiments.
    :rtype: list(int)
    """

    # using the custom API in the sequencing_facility app
    response = uploader.do_get_request(
            '%s_experiment' % uploader.tardis_app_name,
            {'schema_namespace': schema_namespace,
             'parameter_name': 'run_id',
             # 'parameter_type': 'string',
             'parameter_value': run_id})
    if response.ok:
        data = response.json()
        if 'objects' in data:
            return [o.get('id', None) for o in data['objects']]
    else:
        raise response.raise_for_status()


# TODO: Once MyTardis develop supports it, we can use the
#       uploader.query_objectacl method instead of this
def query_objectacl(uploader, object_id, content_type='experiment',
                    acl_ownership_type=u'Owner-owned'):

    acl_ownership_type = uploader._get_ownership_int(acl_ownership_type)

    query_params = {u'object_id': object_id,
                    u'content_type': content_type,
                    u'aclOwnershipType': acl_ownership_type}

    response = uploader.do_get_request(
            '%s_objectacl' % uploader.tardis_app_name,
            query_params)
    return response.json()


def trash_experiments_server(uploader, experiment_ids):
    """

    :param uploader: Uploader object instance.
    :type uploader: MyTardisUploader
    :param experiment_ids: A list of experiement IDs
    :type experiment_ids: list(int)
    :return:
    :rtype:
    """

    for expt in experiment_ids:
        url_template = urljoin(uploader.mytardis_url,
                               '/apps/' + uploader.tardis_app_name + '/api/%s')

        response = uploader._do_request('PUT', 'trash_experiment/%s' % expt,
                                        api_url_template=url_template)
        if response.ok:
            logger.info('Moved experiment %s to trash' % expt)
        else:
            logger.info('Error trying to move experiment %s to trash' % expt)
            raise response.raise_for_status()


def create_experiment_on_server(experiment, uploader):
    """

    :type experiment: mytardis_models.Experiment
    :type uploader: mytardis_uploader.MyTardisUploader
    :return: The url path of the experiment created.
    :rtype: str
    """
    expt_url = uploader.create_experiment(experiment.package())

    return expt_url


def format_instrument_name(instrument_model, instrument_id):
    """

    :type instrument_model: str
    :type instrument_id: str
    :rtype: unicode
    """
    instrument_name = u'%s %s' % (instrument_model, instrument_id)
    return instrument_name


def create_fastq_dataset_on_server(
        proj_id, proj_expt, experiments, uploader,
        run_expt_link=None,
        project_expt_link=None,
        fastqc_summary=None):
    """

    :type proj_expt: Experiment
    :type experiments: list[str]
    :type uploader: mytardis_uploader.MyTardisUploader
    :rtype: str
    """
    fastq_dataset = Dataset()
    fastq_dataset.experiments = experiments
    end_date = _format_day(proj_expt.end_time)
    fastq_dataset.description = 'FASTQ reads, %s, %s' % (proj_id, end_date)

    dataset_params = NucleotideRawReadsDataset()
    dataset_params.from_dict(proj_expt.parameters.to_dict(),
                             existing_only=True)
    dataset_params.ingestor_useragent = uploader.user_agent

    if run_expt_link is not None:
        dataset_params.run_experiment = run_expt_link
    if project_expt_link is not None:
        dataset_params.project_experiment = project_expt_link
    fastq_dataset.parameter_sets.append(dataset_params)

    # This second (hidden) parameter_set, provides a summary of
    # FastQC results for every sample in the project, used for
    # rendering and overview table
    if fastqc_summary is not None:
        fqc_version = fastqc_summary.get('fastqc_version', '')
        fastqc_summary_params = HiddenFastqcProjectSummary()
        fastqc_summary_params.hidden_fastqc_summary_json = \
            json.dumps(fastqc_summary)
        fastqc_summary_params.fastqc_version = fqc_version
        fastq_dataset.parameter_sets.append(fastqc_summary_params)

    instrument_name = format_instrument_name(
        proj_expt.parameters.instrument_model,
        proj_expt.parameters.instrument_id)

    return uploader.create_dataset(fastq_dataset.package(),
                                   instrument=instrument_name)


def create_run_config_dataset_on_server(run_expt, run_expt_url, uploader):
    run_id = run_expt.parameters.run_id
    config_dataset = Dataset()
    config_dataset.experiments = [run_expt_url]
    config_dataset.description = 'Configuration and logs for %s' % run_id
    config_params = IlluminaRunConfig()
    config_params.run_id = run_id
    config_dataset.parameters = config_params
    instrument_name = format_instrument_name(
        run_expt.parameters.instrument_model,
        run_expt.parameters.instrument_id)
    return uploader.create_dataset(config_dataset.package(),
                                   instrument=instrument_name)


def register_project_fastq_datafiles(run_id,
                                     proj_path,
                                     samplesheet,
                                     dataset_url,
                                     uploader,
                                     fastqc_data=None,
                                     fast_mode=False):

    sample_dict = samplesheet_to_dict_by_id(samplesheet)

    # Upload datafiles for the FASTQ reads in the project, for each Sample
    for sample_path, sample_name in get_sample_directories(proj_path):
        for fastq_path in get_fastq_read_files(sample_path):

                sample_id = get_sample_id_from_fastq_filename(fastq_path)

                # sample_name may not be in the SampleSheet dict
                # if we are dealing the unaligned reads (eg sample_names
                # lane1, lane2 etc)
                # So, we grab either the info for the sample_name, or
                # we use an empty dict (where all params will then be
                # the default missing values)

                # TODO: ensure we can also deal with unbarcoded runs,
                #       where Index is "NoIndex" and sample_names are
                #       lane1, lane2 etc.

                sampleinfo = sample_dict.get(sample_name, {})

                reference_genome = sampleinfo.get('SampleRef', '')
                index_sequence = sampleinfo.get('Index', '')
                is_control = sampleinfo.get('Control', '')
                recipe = sampleinfo.get('Recipe', '')
                operator = sampleinfo.get('Operator', '')
                description = sampleinfo.get('Description', '')
                project = sampleinfo.get('SampleProject', '')
                lane = sampleinfo.get('Lane', None)
                if lane is not None:
                    lane = int(lane)
                # flowcell ID is already attached to the Dataset metadata
                # and can be inferrd from the sample_id, so we don't make
                # a special field for it
                # fcid = sampleinfo.get('FCID')

                # the read number isn't encoded in SampleSheet.csv, so we
                # extract it from the FASTQ filename instead
                read = parse_fastq_filename(fastq_path).get('read', None)

                parameters = {'run_id': run_id,
                              'sample_id': sample_id,
                              'sample_name': sample_name,
                              'reference_genome': reference_genome,
                              'index_sequence': index_sequence,
                              'is_control': is_control,
                              'recipe': recipe,
                              'operator_name': operator,
                              'description': description,
                              'project': project,
                              'lane': lane,
                              'read': read,
                              }

                fqc_completed_list = []
                if fastqc_data is not None:
                    fqc_completed_list = [s['filename']
                                          for s in fastqc_data['samples']]

                filename = os.path.basename(fastq_path)
                if filename in fqc_completed_list:
                    # grab the FastQC data for just this FASTQ file
                    sample_fqcdata = next((s for s in fastqc_data['samples'] if
                                           s['filename'] == filename), None)
                    basic_stats = sample_fqcdata['basic_stats']
                    parameters.update(basic_stats)
                elif not fast_mode:
                    # If there is no FastQC data with read counts etc for
                    # this sample (e for Undetermined_indicies) we calculate
                    # our own
                    logger.info("Calculating number of reads for: %s",
                                fastq_path)
                    parameters['number_of_reads'] = \
                        get_number_of_reads_fastq(fastq_path)
                    logger.info("Calculating read length for: %s", fastq_path)
                    parameters['read_length'] = \
                        get_read_length_fastq(fastq_path)

                fq_datafile = DataFile()
                datafile_params = FastqRawReads()
                datafile_params.from_dict(parameters, existing_only=True)
                fq_datafile.parameters = datafile_params
                datafile_parameter_sets = fq_datafile.package_parameter_sets()

                replica_url = fastq_path
                if uploader.storage_mode == 'shared':
                    replica_url = get_shared_storage_replica_url(
                        uploader.storage_box_location,
                        fastq_path)

                if fast_mode:
                    md5_checksum = ''  # '__undetermined__'
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
                except Exception, ex:
                    logger.error("Failed to register Datafile: "
                                 "%s", fastq_path)
                    logger.debug("Exception: %s", ex)
                    raise ex

                logger.info("Added Datafile: %s (%s)",
                            fastq_path,
                            dataset_url)


def parse_fastq_filename(filepath):
    """
    Takes a string like:
    /some/path/15-05065_ACCTCA_L001_R1_001.fastq.gz

    Returns a dictionary like:
    {'read': 1,
    'index': u'ACCTCA',
    'set_number': 1,
    'lane': 1,
    'sample_name': u'15-05065'}

    :type filepath: str
    :rtype: dict
    """
    noext = os.path.basename(filepath).split('.fastq.gz')[0]
    sample_name, index, lane, read, set_number = noext.rsplit('_', 4)
    return dict(sample_name=sample_name,
                index=index,
                lane=int(lane[1:]),
                read=int(read[1:]),
                set_number=int(set_number))


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
    parts = parse_fastq_filename(filepath)
    return parts['sample_name']


def get_sample_id_from_fastqc_zip_filename(filepath):
    return os.path.basename(filepath).split('_fastqc.zip')[0]


def get_sample_name_from_fastqc_filename(fastqc_zip_path):
    return os.path.basename(fastqc_zip_path).split('_')[:1]


def generate_fastqc_report_filename(sample_id):
    return sample_id + "_fastqc.html"


def register_project_fastqc_datafiles(run_id,
                                      proj_id,
                                      fastqc_out_dir,
                                      dataset_url,
                                      uploader,
                                      fast_mode=False):

    # Upload datafiles for the FASTQC output files
    for fastqc_zip_path in get_fastqc_zip_files(fastqc_out_dir):

            fastqc_version = \
                parse_fastqc_data_txt(fastqc_zip_path)['fastqc_version']
            sample_id = get_sample_name_from_fastqc_filename(fastqc_zip_path)
            parameters = {'run_id': run_id,
                          'project': proj_id,
                          'sample_id': sample_id,
                          'fastqc_version': fastqc_version}

            fqc_datafile = DataFile()
            fqc_params = FastqcOutput()
            fqc_params.from_dict(parameters)
            fqc_datafile.parameters = fqc_params
            datafile_parameter_sets = fqc_datafile.package_parameter_sets()

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
            except Exception, ex:
                logger.error("Failed to register Datafile: "
                             "%s", fastqc_zip_path)
                logger.debug("Exception: %s", ex)
                raise ex

            logger.info("Added Datafile: %s (%s)",
                        fastqc_zip_path,
                        dataset_url)


def upload_fastqc_reports(fastqc_out_dir, dataset_url, uploader):
    """
    Uploads the HTML reports generated by FastQC, after extracting them
    from the zipped output produced by the program. Converts reports to be
    self-contained with inline Base64 encoded images if it's not already the
    case (as with older FastQC versions).

    We need to host 'live' html reports at a location that will always be
    immediately available. For this reason, the storage box associated with
    the provided uploader should be writable, and a location that doesn't get
    archived to 'offline' tape storage with delayed access.

    :param fastqc_out_dir:
    :type fastqc_out_dir:
    :param dataset_url:
    :type dataset_url:
    :param uploader:
    :type uploader:
    :return:
    :rtype:
    """

    for fastqc_zip_path in get_fastqc_zip_files(fastqc_out_dir):
        fastqc_data_tables = parse_fastqc_data_txt(fastqc_zip_path)
        fqc_version = fastqc_data_tables['fastqc_version']

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
            sample_id = get_sample_id_from_fastqc_zip_filename(fastqc_zip_path)
            report_dir = join(tmp_path, sample_id + "_fastqc")
            report_file = join(report_dir, 'fastqc_report.html')
            inline_report_filename = generate_fastqc_report_filename(sample_id)
            inline_report_abspath = join(report_dir, inline_report_filename)
            if LooseVersion(fqc_version) < LooseVersion('0.11.3'):
                # convert fastqc_report.html to version with inline images
                from standalone_html import make_html_images_inline
                make_html_images_inline(report_file, inline_report_abspath)
            else:
                os.rename(report_file, inline_report_abspath)

            uploader.upload_file(inline_report_abspath, dataset_url)
            logger.info("Added Datafile (FastQC report): %s (%s)",
                        inline_report_abspath,
                        dataset_url)


def get_fastqc_summary_for_project(fastqc_out_dir):
    """

    Returns a nested data structure representing FastQC results
    for all samples in a project. Includes QC PASS/WARN/FAIL tests,
    original FASTQ filename and some basic statistics (eg number of reads,
    GC content).

    eg:

    {"samples" : [  { "sample_id": "OCT4-15_TGGTGA_L001_R1_001",
                      "sample_name": "OCT4-15",
                      "index": "TGGTGA",
                      "lane": "1",
                      "read_pair": "1",
                      "illumina_sample_sheet": { .... },
                      "qc_checks" : [("Basic Statistics", "PASS"),],
                      "basic_stats": {"number_of_reads": 1042034623, ...},
                      "filename" : "OCT4-15_TGGTGA_L001_R1_001.fastq.gz",
                      "fastqc_report_filename":
                                   "OCT4-15_TGGTGA_L001_R1_001_fastqc.html",
                      ...
                    }.
                    { "sample_id": .... },
                   ...
                 ],
    "fastqc_version" : "0.11.2",
    }

    :type fastqc_out_dir: str
    :rtype project_summary: dict
    """
    project_summary = {u'samples': [], u'fastqc_version': None}
    fqc_detailed_data = {}
    for fastqc_zip_path in get_fastqc_zip_files(fastqc_out_dir):
        qc_pass_fail_table_raw = parse_fastqc_summary_txt(fastqc_zip_path)
        fastq_filename = qc_pass_fail_table_raw[0][2]
        fqfile_details = parse_fastq_filename(fastq_filename)
        sample_id = get_sample_id_from_fastq_filename(fastq_filename)
        fqc_report_filename = generate_fastqc_report_filename(sample_id)

        qc_pass_fail_table = []
        # remove the filename column
        for result, check, fn in qc_pass_fail_table_raw:
            qc_pass_fail_table.append((check, result))

        # project_summary.append((sample_id, qc_pass_fail_table))
        fqc_detailed_data[sample_id] = parse_fastqc_data_txt(fastqc_zip_path)
        basic_stats = _extract_fastqc_basic_stats(fqc_detailed_data, sample_id)
        sample_data = {u'sample_id': sample_id,
                       u'sample_name': fqfile_details['sample_name'],
                       u'qc_checks': qc_pass_fail_table,
                       u'basic_stats': basic_stats,
                       u'filename': fastq_filename,
                       u'fastqc_report_filename': fqc_report_filename,
                       u'index': fqfile_details['index'],
                       u'lane': fqfile_details['lane'],
                       u'read': fqfile_details['read'],
                       # TODO: include the SampleSheet line for this sample
                       u'illumina_sample_sheet': { },
                       }
        project_summary[u'samples'].append(sample_data)
        project_summary[u'fastqc_version'] = \
            fqc_detailed_data[sample_id]['fastqc_version']

    return project_summary

# def get_project_directories(bcl2fastq_out_dir):
#     """
#     Returns an iterator that gives "Project_" directories from
#     the bcl2fastq output directory.
#
#     :type bcl2fastq_out_dir: str
#     :return: Iterator
#     """
#     for item in os.listdir(bcl2fastq_out_dir):
#         proj_path = join(bcl2fastq_out_dir, item)
#         # we look for directories named Project_*
#         # NOTE: an alternative method would be to use the projects explicitly
#         # listed in SameleSheet.csv, treating it as the one source of
#         # truth and avoiding cases where someone decided to create
#         # an additional Project_ directory (eg Project_Bla.old)
#         if isdir(proj_path) and ('Project_' in item):
#             yield proj_path


def get_demultiplexer_info(demultiplexed_output_path):
    """
    Determine which demultiplexing program (eg bcl2fastq) and commandline
    options that were used to partition reads from the run into individual
    samples (based on index).

    eg. {'version': 'bcl2fastq 1.8.3,
         'commandline_options':
         '--input-dir ./Data/Intensities/BaseCalls ' +
         '--output-dir ./130613_SNL177_0029_AH0EPTADXX.pc ' +
         '--sample-sheet ./SampleSheet.csv --no-eamss'}

    :param demultiplexed_output_path: bcl2fastq (or similar) output path
    :type demultiplexed_output_path: str
    :return: Name and version number of program used for demultiplexing reads.
    :rtype: dict
    """
    # Parse DemultiplexConfig.xml for bcl2fastq version
    # NOTE: parsing Makefile would be another option
    demulti_config_path = join(demultiplexed_output_path,
                               "DemultiplexConfig.xml")
    if exists(demulti_config_path):
        with open(demulti_config_path, 'r') as f:
            xml = xmltodict.parse(f)
            version = xml['DemultiplexConfig']['Software']['@Version']
            version = ' '.join(version.split('-'))
            cmdline = xml['DemultiplexConfig']['Software']['@CmdAndArgs']
            cmdline = cmdline.split(' ', 1)[1]
            return {'version': version,
                    'commandline_options': cmdline}

    return {}
    """
    # get the verison assuming the local version of bcl2fastq was used
    out = subprocess.check_output("bcl2fastq --version", shell=True)
    if len(out) >= 2 and 'bcl2fastq' in out[1]:
        version = out[1].strip()
        return {'version': version,
                'commandline_options': None}

    # get the version of locally installed tagdust
    out = subprocess.check_output("tagdust --version", shell=True)
    if len(out) >= 1 and 'Tagdust' in out[1]:
        version = out[1].strip()
        return {'version': version,
                'commandline_options': None}
    """


def get_bcl2fastq_output_dir(options, run_id, run_path):
    return options.bcl2fastq_output_path.format(
        run_id=run_id,
        run_path=run_path
    )


def get_run_id_from_path(run_path):
    return os.path.basename(run_path.strip(os.path.sep))


def get_sample_directories(project_path):
    """
    Returns an iterator that gives tuples of ("Sample_" directory, Sample ID)
    from within a bcl2fastq Project_ directory.

    :type project_path: str
    :return: Iterator
    """
    for item in os.listdir(project_path):
        sample_path = join(project_path, item)
        if isdir(sample_path) and ('Sample_' in item):
            yield sample_path, item.split('Sample_')[1]


def get_fastq_read_files(sample_path):
    """
    Returns an iterator that gives gzipped FASTQ read files from a
    a bcl2fastq Project_*/Sample_* directory.

    :type sample_path: str
    :return: Iterator
    """
    for item in os.listdir(sample_path):
        fastq_path = join(sample_path, item)

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
           isfile(fastq_path):
            yield fastq_path


def get_fastqc_zip_files(fastqc_out_path):
    """
    Returns an iterator that gives zipped FASTQC result files.

    :type fastqc_out_path: str
    :return: Iterator
    """
    for item in os.listdir(fastqc_out_path):
        fastqc_path = join(fastqc_out_path, item)

        # we just return things with the .zip extension
        if os.path.splitext(item)[-1:][0] == '.zip' and \
           isfile(fastqc_path):
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


def run_fastqc_on_project(proj_path,
                          output_directory=None,
                          fastqc_bin=None,
                          threads=2):

    if output_directory is None:
        output_directory = get_fastqc_output_directory(proj_path)
        if exists(output_directory):
            logger.error("FastQC - output directory already exists: %s",
                         output_directory)
            return None
        try:
            os.mkdir(output_directory)
        except OSError:
            logger.error("FastQC - couldn't create output directory: %s",
                         output_directory)
            return None

    if (not exists(output_directory)) or \
       (not isdir(output_directory)):
        logger.error("FastQC - output path %s isn't a directory "
                     "or doesn't exist.",
                     output_directory)
        return None

    fastq_files = []
    for sample_path, sample_id in get_sample_directories(proj_path):
        fastq_files.extend([f for f in get_fastq_read_files(sample_path)])

    fqc_output_directory = run_fastqc(fastq_files,
                                      output_directory=output_directory,
                                      fastqc_bin=fastqc_bin,
                                      threads=threads)
    return fqc_output_directory


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
    def _parse(fh):
        summary = [tuple(line.strip().split('\t')) for line in fh]
        return summary

    return parse_file_from_zip(zip_file_path,
                               'summary.txt',
                               _parse)


def parse_fastqc_data_txt(zip_file_path):
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

    return parse_file_from_zip(zip_file_path,
                               'fastqc_data.txt',
                               _parse)


def _extract_fastqc_basic_stats(fastqc_data, sample_id):
    basic_stats = fastqc_data[sample_id]['Basic Statistics']
    stats = {}
    for row in basic_stats['rows']:
        k = row[0]
        v = row[1]
        if 'Total Sequences' in k:
            stats['number_of_reads'] = int(v)
        if 'Sequences flagged as poor quality' in k:
            stats['number_of_poor_quality_reads'] = int(v)
        if 'Sequence length' in k:
            stats['read_length'] = int(v)
        if '%GC' in k:
            stats['percent_gc'] = float(v)

    return stats


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

    :type storage_box_location: str
    :type file_path: str
    :rtype: str
    """
    if storage_box_location:
        replica_url = os.path.relpath(os.path.normpath(file_path),
                                      storage_box_location)
        return replica_url


def create_tmp_dir(*args, **kwargs):
    tmpdir = mkdtemp(*args, **kwargs)
    global TMPDIRS
    TMPDIRS.append(tmpdir)
    return tmpdir


def get_mytardis_seqfac_app_version(uploader):
    """
    Make a REST call to /apps/sequencing-facility/version to
    determine the version of the sequencing-facility app being
    used by MyTardis.

    :param uploader: The uploader instance
    :type uploader: :py:class:`mytardis_uploader.MyTardisUploader`
    :return: The version string of the remote app
    :rtype: str
    """
    url_template = urljoin(uploader.mytardis_url,
                           '/apps/' + uploader.tardis_app_name + '/api/%s')
    response = uploader._do_request('GET', 'version',
                                    api_url_template=url_template)
    d = response.json()
    version = d.get('version', None)
    return version


def dump_schema_fixtures_as_json():
    fixtures = []
    for name, klass in models.__dict__.items():
        if (not name.startswith('__') and
                issubclass(klass, models.MyTardisParameterSet)):
            fixtures.extend(klass().to_schema())
            fixtures.extend(klass().to_parameter_schema())

    print json.dumps(fixtures, indent=2)


def pre_ingest_checks(options):
    """
    Runs a series of pre-ingestion checks to ensure prerequisite files
    exist, and that the run is ready to ingest and generally sane.

    :param options: An argparse-like options object.
    :type options: :py:class:`object`
    :return: True if all checks pass, False otherwise
    :rtype: bool
    """

    def get_command_stdout(cmd):
        """
        Runs a shell command using subprocess.check_output and returns stdout,
        irrespective of any non-zero exit codes.

        :param cmd: The shell command
        :type cmd: basestring
        :return: stdout generated by the command
        :rtype: basestring
        """
        try:
            stdout = subprocess.check_output(cmd, shell=True)
        except subprocess.CalledProcessError as e:
            stdout = str(e.output)
        return stdout

    run_path = options.path
    run_id = get_run_id_from_path(run_path)

    bcl2fastq_output_dir = get_bcl2fastq_output_dir(options,
                                                    run_id,
                                                    run_path)

    if not exists(bcl2fastq_output_dir):
        logger.error("Aborting - bcl2fastq output directory (%s) not found. "
                     "Is your setting for --bcl2fastq-output-path correct ?" %
                     bcl2fastq_output_dir)
        return False

    if not exists(join(run_path, 'RTAComplete.txt')):
        logger.error("Aborting - 'RTAComplete.txt' not found. "
                     "Is the run still in progress ?")
        return False

    if not exists(join(run_path, 'RunInfo.xml')):
        logger.error("Aborting - 'RunInfo.xml' not found.")
        return False

    config_filename = None
    for filename in os.listdir(join(run_path, 'Config')):
        if "_Effective.cfg" in filename:
            config_filename = filename
    if not config_filename:
        logger.error("Aborting - cannot find Config/*_Effective.cfg file.")
        return False

    if not exists(join(run_path, 'SampleSheet.csv')):
        logger.error("Aborting - 'SampleSheet.csv' not found.")
        return False

    try:
        samplesheet, chemistry = parse_samplesheet(join(run_path,
                                                        'SampleSheet.csv'))
    except IOError:
        logger.error("Aborting - unable to parse SampleSheet.csv file.")
        return False

    projects = get_project_ids_from_samplesheet(
        samplesheet,
        include_no_index_name='Undetermined_indices')

    for p in projects:
        p_path = join(bcl2fastq_output_dir, 'Project_%s' % p)
        if p == 'Undetermined_indices':
            p_path = join(bcl2fastq_output_dir, 'Undetermined_indices')

        if not exists(p_path):
            logger.error("Aborting - project directory '%s' is missing.",
                         p_path)
            return False
            # sys.exit(1)

    # Check for zero sized .bcl files, indicative of a failed RTA transfer
    if options.ignore_zero_sized_bcl_check:
        zero_size_files = None
    else:
        bcl_path = join(run_path, 'Data/Intensities/BaseCalls')
        zero_size_files = get_command_stdout(
            'find %s -size 0 -print | grep ".bcl$"' % bcl_path)

    if zero_size_files.strip():
        logger.error("Aborting - some BCL files in %s are zero size."
                     "As a result, there may be issues with the demultiplexed "
                     "FASTQ files. Use the flag --ignore-zero-sized-bcl-check "
                     "if you really want to proceed.")

    # Check for keyword '[fF]ail' in <run_path>/Data/RTALogs/*, indicative
    # of a failed transfer
    fail_logs = get_command_stdout(
        'grep [fF]ail %s/*' % join(run_path, 'Data/RTALogs'))

    if fail_logs.strip():
        logger.warn("WARNING - logs in Data/RTALogs contain failure messages."
                    "There may have been issues transfering the "
                    "basecalling data (Data/Intensities/BaseCalls/*/*.bcl) from"
                    "the instrument.",
                     p_path)
        return False

    return True


def ingest_run(run_path=None):
    global logger
    logger = setup_logging()
    global TMPDIRS
    TMPDIRS = []

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
        argparser.add_argument('--dump-fixtures',
                               dest='dump_fixtures',
                               action='store_true')
        argparser.add_argument('--live-storage-box-name',
                               dest='live_storage_box_name',
                               default='live',
                               type=str,
                               metavar='LIVE_STORAGE_BOX_NAME')
        argparser.add_argument('--replace-duplicate-runs',
                               dest='replace_duplicate_runs',
                               type=bool,
                               default=False,
                               metavar='REPLACE_DUPLICATE_RUNS')
        argparser.add_argument('--trash-user',
                               dest='trash_user',
                               default='__trashman__',
                               type=str,
                               metavar='TRASH_USER')
        argparser.add_argument('--trash-group',
                               dest='trash_group',
                               default='__trashcan__',
                               type=str,
                               metavar='TRASH_GROUP')
        argparser.add_argument('--ignore-zero-sized-bcl-check',
                               dest='ignore_zero_sized_bcl_check',
                               type=bool,
                               default=False,
                               metavar='IGNORE_ZERO_SIZED_BCL_CHECK')

    parser, options = get_config(add_extra_options_fn=extra_config_options)

    if options.dump_fixtures:
        dump_schema_fixtures_as_json()
        sys.exit(1)

    validate_config(parser, options)

    # Before creating any records on the server we first check that certain
    # prerequisite files exist, that the run is complete and is generally in a
    # 'sane' state suitable for ingestion.
    logger.info("Running pre-ingestion checks & validation.")
    if not pre_ingest_checks(options):
        sys.exit(1)

    # exclude_patterns = \
    #     get_exclude_patterns_as_regex_list(options.exclude)

    uploader = MyTardisUploader(
        options.url,
        options.username,
        password=options.password,
        api_key=options.api_key,
        storage_mode=options.storage_mode,
        storage_box_location=options.storage_base_path,
        storage_box_name=options.storage_box_name,
        verify_certificate=options.verify_certificate
    )

    # This uploader instance is associated with a MyTardis storage box
    # to which MyTardis has write access, and is intended for files
    # that are served 'live' (eg HTML reports etc), and shouldn't be
    # archived to tape
    writable_storage_uploader = MyTardisUploader(
        options.url,
        options.username,
        password=options.password,
        api_key=options.api_key,
        storage_mode='upload',
        # storage_box_location='/data/cached',
        # storage_box_name='live',
        # storage_box_name='object_store',
        storage_box_name=options.live_storage_box_name,
        verify_certificate=options.verify_certificate,
    )

    # this custom attribute on the uploader is the name of the
    # sequencing_facility MyTardis app, used to contruct the URLs for
    # some app-specific REST API calls
    uploader.tardis_app_name = 'sequencing-facility'

    ingestor_version = SemanticVersion(mytardis_uploader.__version__)
    logger.info("Verifying MyTardis server app '%s' matches the ingestor "
                "version (%s)." % (uploader.tardis_app_name, ingestor_version))
    seqfac_app_version = SemanticVersion(
        get_mytardis_seqfac_app_version(uploader))
    if seqfac_app_version != ingestor_version:
        logger.error("Ingestor (%s) / server (%s) version mismatch." %
                     (seqfac_app_version, mytardis_uploader.__version__))
        raise Exception("Version mismatch.")

    # Create an Experiment representing the overall sequencing run
    if not run_path:
        run_path = options.path
    run_expt = create_run_experiment_object(run_path)
    # run_id = get_run_id_from_path(run_path)
    run_id = run_expt.parameters.run_id

    duplicate_runs = get_experiments_from_server_by_run_id(
        uploader, run_id,
        'http://www.tardis.edu.au/schemas/ngs/run/illumina')
    duplicate_projects = get_experiments_from_server_by_run_id(
        uploader, run_id,
        'http://www.tardis.edu.au/schemas/ngs/project')

    if duplicate_runs or duplicate_projects:
        matching = [unicode(i) for i in duplicate_runs]
        matching += [unicode(i) for i in duplicate_projects]
        logger.warn("Duplicate runs/projects already exist on server: %s (%s)",
                    run_id, ', '.join(matching))
        if options.replace_duplicate_runs:
            trash_experiments_server(uploader, matching)
        else:
            logger.error("Please manually remove existing run before "
                         "ingesting: %s (%s)", run_id, ', '.join(matching))
            raise Exception()

    # The directory where bcl2fastq puts its output,
    # in Project_* directories
    bcl2fastq_output_dir = get_bcl2fastq_output_dir(options,
                                                    run_id,
                                                    run_path)

    run_expt.institution_name = options.institute
    run_expt.description = options.description
    run_expt.parameters.ingestor_useragent = uploader.user_agent
    demultiplexer_info = get_demultiplexer_info(bcl2fastq_output_dir)
    run_expt.parameters.demultiplexing_program = \
        demultiplexer_info.get('version', '')
    run_expt.parameters.demultiplexing_commandline_options = \
        demultiplexer_info.get('commandline_options', '')
    if not run_expt.description:
        run_expt.description = "%s run %s, ingested on %s" % \
                               (run_expt.parameters.instrument_model,
                                run_id,
                                datetime.now().isoformat(' '))

    try:
        run_expt_url = create_experiment_on_server(run_expt, uploader)

        # Take just the path of the experiment, eg /api/v1/experiment/187/
        run_expt_url = urlparse(run_expt_url).path

        for group in options.experiment_owner_groups:
            uploader.share_experiment_with_group(run_expt_url, group)

    except Exception, e:
        logger.error("Failed to create Experiment for sequencing run: %s",
                     run_path)
        logger.error("Exception: %s: %s", type(e).__name__, e)
        raise e

    logger.info("Created Run Experiment: %s (%s)",
                run_id,
                run_expt_url)

    samplesheet_path = join(run_path, 'SampleSheet.csv')
    samplesheet, chemistry = parse_samplesheet(samplesheet_path)

    # Under the Run Experiment we create a Dataset with the IlluminaRunConfig
    # schema containing the SampleSheet.csv, maybe also some logs and
    # config files
    try:
        config_dataset_url = create_run_config_dataset_on_server(run_expt,
                                                                 run_expt_url,
                                                                 uploader)
        config_dataset_url = urlparse(config_dataset_url).path
        uploader.upload_file(samplesheet_path, config_dataset_url)
        logger.info("Created config & logs dataset for sequencing run: %s",
                    config_dataset_url)
    except Exception as e:
        logger.error("Failed to create config & logs dataset for sequencing "
                     "run: %s",
                     run_path)
        logger.error("Exception: %s: %s", type(e).__name__, e)
        raise e

    projects = get_project_ids_from_samplesheet(
        samplesheet,
        include_no_index_name='Undetermined_indices')

    for proj_id in projects:
        proj_expt = create_project_experiment_object(
            proj_id,
            run_expt,
            run_expt_link=run_expt_url)

        proj_expt.parameters.ingestor_useragent = uploader.user_agent
        proj_expt.parameters.demultiplexing_program = \
            run_expt.parameters.demultiplexing_program
        proj_expt.parameters.demultiplexing_commandline_options = \
            run_expt.parameters.demultiplexing_commandline_options

        proj_path = join(bcl2fastq_output_dir, 'Project_' + proj_id)

        if proj_id == 'Undetermined_indices':
            proj_path = join(bcl2fastq_output_dir, proj_id)

        fastqc_out_dir = get_fastqc_output_directory(proj_path)

        # Run FastQC if output doesn't exist.
        # We output to a tmp dir since we don't expect to have write
        # permissions to the primary data directory
        if options.run_fastqc \
                and not options.fast \
                and not exists(fastqc_out_dir):
            fqc_tmp_dir = create_tmp_dir()
            fastqc_out_dir = run_fastqc_on_project(
                proj_path,
                output_directory=fqc_tmp_dir,
                fastqc_bin=options.fastqc_bin,
                threads=int(options.threads)
            )

        fqc_summary = {}
        if exists(fastqc_out_dir):
            fqc_summary = get_fastqc_summary_for_project(fastqc_out_dir)

        # Create an Experiment for each real Project in the run
        # (except Undetermined_indices Datasets which will only be associated
        #  with the parent 'run' Experiment)
        project_url = None
        if proj_id != 'Undetermined_indices':
            try:
                project_url = create_experiment_on_server(proj_expt, uploader)

                project_url = urlparse(project_url).path

                for group in options.experiment_owner_groups:
                    uploader.share_experiment_with_group(project_url, group)

            except Exception, e:
                logger.error("Failed to create Experiment for project: %s",
                             proj_id)
                logger.debug("Exception: %s", e)
                raise e

            logger.info("Created Project Experiment: %s (%s)",
                        project_url,
                        proj_id)

        # samples_in_project = [s['SampleID'] for s in samplesheet]
        # sample_desc = 'FASTQ reads for samples '+', '.join(samples_in_project)

        # We associate Undetermined_indices Datasets with the overall 'run'
        # Experiment only (unlike proper Project Datasets which also have
        # their own Project Experiment).
        if proj_id == 'Undetermined_indices':
            parent_expt_urls = [run_expt_url]
        else:
            parent_expt_urls = [project_url, run_expt_url]

        ############################################
        # Create the FastQC Dataset for the project
        if exists(fastqc_out_dir):
            try:
                # TODO: fq_dataset_url should actually be a URL .. but ..
                # we have a chicken-egg problem here - we want the URL
                # for the FASTQ dataset here, to add as a parameter, but
                # we are adding the FastQC dataset URL so we can add it's
                # URL to the FASTQ dataset.
                # We need to be able to update the FastQC dataset parameters
                # in a second API call after we've added the FASTQ dataset.
                # fq_dataset_url = "%s__%s" % (run_id, proj_id)
                fq_dataset_url = project_url # placeholder

                # Then discard parts, repopulate some parameters
                fqc_dataset = create_fastqc_dataset_object(
                    run_id,
                    proj_id,
                    proj_expt.end_time,
                    parent_expt_urls,
                    fq_dataset_url,
                    fastqc_summary=fqc_summary)

                fqc_dataset.parameters.ingestor_useragent = uploader.user_agent

                # TODO: encapsulate this in a create_fastqc_dataset function
                #       for consistency ?
                fqc_dataset_url = uploader.create_dataset(fqc_dataset.package())

                # Take just the path, eg: /api/v1/dataset/363
                fqc_dataset_url = urlparse(fqc_dataset_url).path
                # Add the LINK parameter from the FASTQ dataset to it's
                # associated FastQC dataset
                proj_expt.parameters.fastqc_dataset = fqc_dataset_url

            except Exception, e:
                logger.error("Failed to create FastQC Dataset for Project: %s",
                             proj_id)
                logger.debug("Exception: %s", e)
                raise e

            logger.info("Created FastQC Dataset: %s (%s)",
                        fqc_dataset_url,
                        proj_id)

            # We don't add the FastQC zips for those in temporary directories
            # eg, for 'Undetermined_indices' (since these won't be present on
            # shared storage).
            # We always upload the html reports to be serverd live (below).
            if fastqc_out_dir and fastqc_out_dir not in TMPDIRS:
                register_project_fastqc_datafiles(run_id,
                                                  proj_id,
                                                  fastqc_out_dir,
                                                  fqc_dataset_url,
                                                  uploader,
                                                  fast_mode=options.fast)

            upload_fastqc_reports(fastqc_out_dir, fqc_dataset_url,
                                  writable_storage_uploader)

        #################################################################
        # Create the FASTQ Dataset for the project, associated with both
        # the overall run Experiment, and the project Experiment
        try:
            fq_dataset_url = create_fastq_dataset_on_server(
                proj_id,
                proj_expt,
                parent_expt_urls,
                uploader,
                fastqc_summary=fqc_summary)
        except Exception, e:
            logger.error("Failed to create Dataset for Project: %s",
                         proj_id)
            logger.debug("Exception: %s", e)
            raise e

        # Take just the path, eg: /api/v1/dataset/363
        fq_dataset_url = urlparse(fq_dataset_url).path

        logger.info("Created FASTQ Dataset: %s (%s)", fq_dataset_url, proj_id)

        try:
            # Create a temporary SampleSheet.csv containing only lines for the
            # current Project, to be uploaded to the FASTQ Dataset
            tmp_dir = create_tmp_dir()
            project_samplesheet_path = join(tmp_dir, 'SampleSheet.csv')
            with open(project_samplesheet_path, 'w') as f:
                lines = filter_samplesheet_by_project(samplesheet_path,
                                                      proj_id)
                f.writelines(lines)

            writable_storage_uploader.upload_file(project_samplesheet_path,
                                                  fq_dataset_url)
            if exists(tmp_dir):
                shutil.rmtree(tmp_dir)

            logger.info("Uploaded SampleSheet.csv for Project: %s (%s)",
                        fq_dataset_url,
                        proj_id)
        except Exception as e:
            logger.error("Uploading SampleSheet.csv for Project failed: %s (%s)",
                        fq_dataset_url,
                        proj_id)

        register_project_fastq_datafiles(run_id,
                                         proj_path,
                                         samplesheet,
                                         fq_dataset_url,
                                         uploader,
                                         fastqc_data=fqc_summary,
                                         fast_mode=options.fast)

    logger.info("Ingestion of run %s complete !", run_id)


@atexit.register
def _cleanup_tmp():
    global TMPDIRS
    for tmpdir in TMPDIRS:
        if exists(tmpdir):
            shutil.rmtree(tmpdir)
            logger.info("Removed temp directory: %s", tmpdir)

if __name__ == "__main__":
    MyTardisUploader.user_agent_name = os.path.basename(sys.argv[0])
    try:
        ingest_run()
    except Exception, e:
        if e != SystemExit:
            import traceback
            # traceback.print_exc(file=sys.stdout)
            logger.debug((traceback.format_exc()))
        _cleanup_tmp()
        logger.error("Ingestion failed.")
        sys.exit(1)

    # since atexit doesn't seem to work
    _cleanup_tmp()

    sys.exit(0)
