#!/usr/bin/env python
from __future__ import print_function, absolute_import, division

__author__ = 'Andrew Perry <Andrew.Perry@monash.edu.au>'

from builtins import (bytes, str, open, super, range,
                      zip, round, input, int, pow, object)
import six
from six.moves.urllib.parse import urlparse, urljoin
import sys
import shutil
from datetime import datetime
import subprocess
from tempfile import mkdtemp
import atexit

import os
from os.path import join, splitext, exists, isdir, isfile
import json

from semantic_version import Version as SemanticVersion
from distutils.version import LooseVersion

from standalone_html import make_html_images_inline

import mytardis_uploader
from mytardis_uploader import MyTardisUploader
from mytardis_uploader import setup_logging, get_config, validate_config
# from mytardis_ngs_ingestor import get_exclude_patterns_as_regex_list

import models
from mytardis_models import Experiment, Dataset, DataFile
from models import DemultiplexedSamplesBase, FastqcOutputBase, \
    FastqcReportsBase, FastqRawReadsBase, HiddenFastqcProjectSummaryBase, \
    IlluminaSequencingRunBase, NucleotideRawReadsDatasetBase, \
    IlluminaRunConfigBase

from illumina.run_info import parse_samplesheet, samplesheet_to_dict, \
    get_project_ids_from_samplesheet, get_number_of_reads_fastq, \
    get_read_length_fastq, rta_complete_parser, runinfo_parser, \
    illumina_config_parser, get_run_id_from_path, get_demultiplexer_info, \
    get_sample_directories, get_fastq_read_files, \
    get_sample_id_from_fastq_filename, get_sample_name_from_fastq_filename, \
    parse_sample_info_from_filename, filter_samplesheet_by_project

# a module level list of temporary directories that have been
# created, so these can be cleaned up upon premature exit
TMPDIRS = []

# we map names from *Effective.cfg to more recognisable names
# (MiSeq seems to be called 'Shock' some places internally)
INSTRUMENT_TYPE_NAMES = {
    'Shock': 'MiSeq'
}

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
    raw_model_name = instrument_config.get('system:instrumenttype', '')
    instrument_model = INSTRUMENT_TYPE_NAMES.get(raw_model_name, raw_model_name)
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
    expt.title = "%s run %s" % (run.instrument_model, run.run_id)
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
    if proj_id:
        project_title = 'Sequencing Project, %s, %s' % (proj_id, end_date)
    else:
        project_title = 'Sequencing Project, %s' % end_date

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

    if not proj_id:
        fqc_dataset_title = "FastQC reports, %s" % end_date

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


def create_fastqc_dataset_on_server(fastqc_dataset, uploader):
    return uploader.create_dataset(fastqc_dataset.package())


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
    if isinstance(schema_namespace, six.string_types):
        query_params[u'name__schema__namespace'] = schema_namespace
    if isinstance(parameter_name, six.string_types):
        query_params[u'name__name'] = parameter_name
    if isinstance(value, (float, int)):
        query_params[u'numeric_value'] = value
    if isinstance(value, six.string_types):
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
    if proj_id:
        fastq_dataset.description = 'FASTQ reads, %s, %s' % (proj_id, end_date)
    else:
        fastq_dataset.description = 'FASTQ reads, %s' % end_date

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
                                     no_sample_directories=False,
                                     fast_mode=False):

    sample_dict = samplesheet_to_dict(samplesheet)

    if no_sample_directories:  # eg, Undetermined_indices
        sample_path_and_name = [(proj_path, None)]
    else:
        sample_path_and_name = get_sample_directories(proj_path)

    # Upload datafiles for the FASTQ reads in the project,
    # for each Sample_ directory
    for sample_path, sample_name in sample_path_and_name:
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

                # the read number isn't encoded in SampleSheet.csv, so we
                # extract it from the FASTQ filename instead
                read = None
                info_from_fn = parse_sample_info_from_filename(fastq_path)
                if info_from_fn is not None:
                    read = info_from_fn.get('read', None)
                    if sample_name is None:
                        sample_name = info_from_fn.get('sample_name', None)
                else:
                    logger.warning("Unrecognized FASTQ filename pattern - "
                                   "skipping: %s", fastq_path)
                    continue

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
                if fastqc_data:
                    fqc_completed_list = [s['filename']
                                          for s in fastqc_data.get('samples',
                                                                   [])
                                          ]

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
                except Exception as ex:
                    logger.error("Failed to register Datafile: "
                                 "%s", fastq_path)
                    logger.debug("Exception: %s", ex)
                    raise ex

                logger.info("Added Datafile: %s (%s)",
                            fastq_path,
                            dataset_url)


def get_sample_id_from_fastqc_zip_filename(filepath):
    return os.path.basename(filepath).split('_fastqc.zip')[0]


def get_sample_name_from_fastqc_filename(fastqc_zip_path):
    return os.path.basename(fastqc_zip_path).split('_')[:1]


def generate_fastqc_report_filename(sample_id):
    return sample_id + "_fastqc.html"


# TODO: Rather than constructing project directories this way
#       we should do generic detection to get a list of project
#       dirs via looking for
#       # v1.8.4
#       {project_id}/Sample_{sample_id}/{sample_id}_{index}_L00{lane}_R{read}_001.fastq.gz or
#       # v2.x
#       {project_id}/{sample_id}/{sample_name}_S{sample_id_num}_L00{lane}_R{read}_001.fastq.gz
#       (where sample_id_num is the index+1 where that sample_id is first seen in the
#        sample sheet)
#
#       eg, if we were to do a recursive search for all *.fastq.gz relative paths
#           take the first part of the path ({project_id}) and make a Set to
#           remove duplicates
def proj_id_to_proj_dir(proj_id, demultiplexer_version_num='1.8.4'):
    if demultiplexer_version_num and \
       LooseVersion(demultiplexer_version_num) >= LooseVersion('2'):
        return proj_id
    else:
        return 'Project_%s' % proj_id


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
            except Exception as ex:
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
                make_html_images_inline(report_file, inline_report_abspath)
            else:
                os.rename(report_file, inline_report_abspath)

            uploader.upload_file(inline_report_abspath, dataset_url)
            logger.info("Added Datafile (FastQC report): %s (%s)",
                        inline_report_abspath,
                        dataset_url)


def get_fastqc_summary_for_project(fastqc_out_dir, samplesheet):
    """

    Returns a nested data structure representing FastQC results
    for all samples in a project. Includes QC PASS/WARN/FAIL tests,
    original FASTQ filename and some basic statistics (eg number of reads,
    GC content). Sorts based on the ordering found in samplesheet
    (SampleSheet.csv).

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

    :param fastqc_out_dir: Path of the output directory containing
                           *_fastqc.zip files.
    :type fastqc_out_dir: str
    :param samplesheet: A list of samplesheet lines as dicts, as output by
                        parse_samplesheet.
    :type samplesheet: list[dict]
    :rtype project_summary: dict
    """
    project_summary = {u'samples': [], u'fastqc_version': None}
    fqc_detailed_data = {}

    def find_sample_index(a):
        """
        Returns the index of a matching sample from SampleSheet.csv, given
        a dictionary of the type returned by parse_sample_info_from_filename.

        :param a: A dictionary representing a Sample ID.
        :type a: dict
        :return: The index of this sample in the SampleSheet.csv, or None if
                 not found
        :rtype: int | None
        """

        # for bcl2fastq2 output files, we have an _S*_ sample number we can use
        if a is None:
            return None
        sample_number = a.get('sample_number', None)
        if sample_number is not None:
            return int(sample_number - 1)

        # otherwise, fall back to finding the sample in the samplesheet
        # to get it's index
        for i, s in enumerate(samplesheet):
            if s['SampleID'] == a['sample_name'] and \
                    (s.get('Index', None) == a['index'] or
                             s.get('index', None) == a['index']) and \
                            int(s['Lane']) == a['lane']:
                return i

        return None

    def sort_fqc_sample_cmp(a, b):
        """
        Sorting comparison (cmp) function to order FastQC output filenames
        to the equivalent ordering in SampleSheet.csv

        :param a: The path to a FastQC output zip file.
        :param b: The path to a FastQC output zip file.
        :type a: str
        :type b: str
        :rtype order: int
        """
        aa = parse_sample_info_from_filename(a, '_fastqc.zip')
        bb = parse_sample_info_from_filename(b, '_fastqc.zip')
        ai = find_sample_index(aa)
        bi = find_sample_index(bb)
        if ai is None or bi is None:
            return 0
        if ai > bi:
            return 1
        if ai < bi:
            return -1
        if ai == bi:
            if aa['lane'] < bb['lane']:
                return 1
            if aa['lane'] > bb['lane']:
                return -1
            if aa['lane'] == bb['lane']:
                if aa['read'] < bb['read']:
                    return 1
                if aa['read'] > bb['read']:
                    return -1
                if aa['read'] == bb['read']:
                    return -1 if aa['set_number'] < bb['set_number'] else 1

        return 0

    fastqc_zips = list(get_fastqc_zip_files(fastqc_out_dir))
    fastqc_zips.sort(cmp=sort_fqc_sample_cmp)

    for fastqc_zip_path in fastqc_zips:
        qc_pass_fail_table_raw = parse_fastqc_summary_txt(fastqc_zip_path)
        fastq_filename = qc_pass_fail_table_raw[0][2]
        fqfile_details = parse_sample_info_from_filename(fastq_filename)
        sample_id = get_sample_id_from_fastq_filename(fastq_filename)
        fqc_report_filename = generate_fastqc_report_filename(sample_id)

        qc_pass_fail_table = []
        # remove the filename column
        for result, check, fn in qc_pass_fail_table_raw:
            qc_pass_fail_table.append((check, result))

        # project_summary.append((sample_id, qc_pass_fail_table))
        fqc_detailed_data[sample_id] = parse_fastqc_data_txt(fastqc_zip_path)
        basic_stats = _extract_fastqc_basic_stats(fqc_detailed_data, sample_id)
        sample_name = fqfile_details.get('sample_name', None)
        lane = fqfile_details.get('lane', None)
        index = fqfile_details.get('index', None)

        # TODO: this could fail if there are two samples of the same name
        #       that share a lane (odd corner case). we should consider
        #       instead using fastq_filename to read the first header of
        #       the actual FASTQ file and extract the index from there
        #       get the index for cases where it isn't in the filename
        if index is None:
            for s in samplesheet:
                if s.get('SampleName', None) == sample_name and \
                   s.get('Lane', None) == str(lane):
                    index = s.get('index', None) or s.get('Index', None)

        sample_data = {u'sample_id': sample_id,
                       u'sample_name': sample_name,
                       u'qc_checks': qc_pass_fail_table,
                       u'basic_stats': basic_stats,
                       u'filename': fastq_filename,
                       u'fastqc_report_filename': fqc_report_filename,
                       u'index': index,
                       u'lane': lane,
                       u'read': fqfile_details.get('read', None),
                       # TODO: include the SampleSheet line for this sample
                       u'illumina_sample_sheet': { },
                       }

        project_summary[u'samples'].append(sample_data)
        project_summary[u'fastqc_version'] = \
            fqc_detailed_data[sample_id]['fastqc_version']

    return project_summary


def get_bcl2fastq_output_dir(options, run_id, run_path):
    return options.bcl2fastq_output_path.format(
        run_id=run_id,
        run_path=run_path
    )


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
    try:
        d = response.json()
    except ValueError as ex:  # JSONDecodeError inherits from ValueError
        logger.error("Invalid response when querying server version.")
        raise ex

    version = d.get('version', None)
    return version


def is_server_version_compatible(ingestor_version, server_version):
    return SemanticVersion(server_version) == SemanticVersion(ingestor_version)


def dump_schema_fixtures_as_json():
    fixtures = []
    for name, klass in models.__dict__.items():
        if (not name.startswith('__') and
                issubclass(klass, models.MyTardisParameterSet)):
            fixtures.extend(klass().to_schema())
            fixtures.extend(klass().to_parameter_schema())

    print(json.dumps(fixtures, indent=2))


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
        if "Effective.cfg" in filename:
            config_filename = filename
    if not config_filename:
        logger.error("Aborting - cannot find Config/*Effective.cfg file.")
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

    if not exists(join(bcl2fastq_output_dir, 'DemultiplexConfig.xml')):
        logger.warning("'DemultiplexConfig.xml' not found.")

    demultiplexer_info = get_demultiplexer_info(bcl2fastq_output_dir)
    demultiplexer_version = demultiplexer_info.get('version', '')

    if not demultiplexer_version:
        logger.warning("Can't determine demultiplexer version")
        # return False

    demultiplexer_version_num = demultiplexer_info.get('version_number', '')
    undetermined_in_root_folder = False
    if demultiplexer_version_num:
        undetermined_in_root_folder = \
            LooseVersion(demultiplexer_version_num) >= LooseVersion('2')

    projects = get_project_ids_from_samplesheet(
        samplesheet,
        include_no_index_name='Undetermined_indices')

    for p in projects:
        p_path = join(bcl2fastq_output_dir,
                      proj_id_to_proj_dir(
                          p,
                          demultiplexer_version_num=demultiplexer_version_num))
        if p == 'Undetermined_indices':
            if undetermined_in_root_folder:
                p_path = bcl2fastq_output_dir
            else:
                p_path = join(bcl2fastq_output_dir, 'Undetermined_indices')

        if not exists(p_path):
            logger.error("Aborting - project directory '%s' is missing.",
                         p_path)
            return False
            # sys.exit(1)

    # Check for zero sized .bcl files, indicative of a failed RTA transfer
    zero_size_files = ''
    if not options.ignore_zero_sized_bcl_check:
        bcl_path = join(run_path, 'Data/Intensities/BaseCalls')
        if exists(bcl_path):
            zero_size_files = get_command_stdout(
                'find %s -size 0 -print | grep ".bcl(.gz)?"' % bcl_path)
        else:
            logger.warning('No Data/Intensities/BaseCalls directory found. '
                           'Skipping zero-sized bcl sanity check.')

    if zero_size_files.strip():
        logger.error("Aborting - some BCL files in %s are zero size."
                     "As a result, there may be issues with the demultiplexed "
                     "FASTQ files. Use the flag --ignore-zero-sized-bcl-check "
                     "if you really want to proceed.")

    # Check for keyword '[fF]ail' in <run_path>/Data/RTALogs/*, indicative
    # of a failed transfer
    rtalogs_path = None
    if exists(join(run_path, 'Data/RTALogs')):
        rtalogs_path = join(run_path, 'Data/RTALogs')
    elif exists(join(run_path, 'RTALogs')):
        rtalogs_path = join(run_path, 'RTALogs')

    if rtalogs_path:
        fail_logs = get_command_stdout(
            'grep [fF]ail %s/*' % rtalogs_path)

        if fail_logs.strip():
            logger.warn("WARNING - logs in Data/RTALogs contain failure "
                        "messages: %s", p_path)
            # return False
    else:
        logger.warn("WARNING - RTALogs or Data/RTALogs directory not found.")

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
                               default='{run_path}/Data/Intensities/BaseCalls',
                               type=str,
                               metavar='BCL2FASTQ_OUTPUT_PATH',
                               help='The path to the bcl2fastq output '
                                    '(fastq.gz files in project/sample '
                                    'directories). The template strings '
                                    '{run_path} and {run_id} can be used to '
                                    'specify a path relative to the run '
                                    'folder, or another path that includes the '
                                    'run_id.')
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
        verify_certificate=options.verify_certificate,
        fast_mode=options.fast,
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
        fast_mode=options.fast,
    )

    # this custom attribute on the uploader is the name of the
    # sequencing_facility MyTardis app, used to contruct the URLs for
    # some app-specific REST API calls
    uploader.tardis_app_name = 'sequencing-facility'

    ingestor_version = mytardis_uploader.__version__
    seqfac_app_version = get_mytardis_seqfac_app_version(uploader)
    logger.info("Verifying MyTardis server app '%s' matches the ingestor "
                "version (%s)." % (uploader.tardis_app_name,
                                   ingestor_version))
    if not is_server_version_compatible(ingestor_version, seqfac_app_version):
        logger.error("Ingestor (%s) / server (%s) version mismatch." %
                     (seqfac_app_version, ingestor_version))
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
        matching = [str(i) for i in duplicate_runs]
        matching += [str(i) for i in duplicate_projects]
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
    demultiplexer_version = demultiplexer_info.get('version', '')
    run_expt.parameters.demultiplexing_program = demultiplexer_version
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

    except Exception as e:
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

    if not demultiplexer_info.get('version', None):
        logger.error("Can't determine demultiplexer version - aborting")
        raise Exception()

    demultiplexer_version_num = demultiplexer_info.get('version_number', '')
    undetermined_in_root_folder = False
    if demultiplexer_version_num:
        undetermined_in_root_folder = \
            LooseVersion(demultiplexer_version_num) >= LooseVersion('2')

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

        proj_path = join(bcl2fastq_output_dir,
                         proj_id_to_proj_dir(
                             proj_id,
                             demultiplexer_version_num=demultiplexer_version_num))

        if proj_id == 'Undetermined_indices':
            if undetermined_in_root_folder:
                proj_path = bcl2fastq_output_dir
            else:
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
        if fastqc_out_dir is not None and exists(fastqc_out_dir):
            fqc_summary = get_fastqc_summary_for_project(fastqc_out_dir,
                                                         samplesheet)

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

            except Exception as e:
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
        if fastqc_out_dir is not None and exists(fastqc_out_dir):
            try:
                # TODO: fq_dataset_url should actually be a URL .. but ..
                # we have a chicken-egg problem here - we want the URL
                # for the FASTQ dataset here, to add as a parameter, but
                # we are adding the FastQC dataset URL so we can add it's
                # URL to the FASTQ dataset.
                # We need to be able to update the FastQC dataset parameters
                # in a second API call after we've added the FASTQ dataset.
                # fq_dataset_url = "%s__%s" % (run_id, proj_id)
                fq_dataset_url = project_url  # placeholder

                # Then discard parts, repopulate some parameters
                fqc_dataset = create_fastqc_dataset_object(
                    run_id,
                    proj_id,
                    proj_expt.end_time,
                    parent_expt_urls,
                    fq_dataset_url,
                    fastqc_summary=fqc_summary)

                fqc_dataset.parameters.ingestor_useragent = uploader.user_agent

                fqc_dataset_url = create_fastqc_dataset_on_server(fqc_dataset,
                                                                  uploader)

                # Take just the path, eg: /api/v1/dataset/363
                fqc_dataset_url = urlparse(fqc_dataset_url).path
                # Add the LINK parameter from the FASTQ dataset to it's
                # associated FastQC dataset
                proj_expt.parameters.fastqc_dataset = fqc_dataset_url

            except Exception as e:
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
        except Exception as e:
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

        if proj_id == 'Undetermined_indices':
            no_sample_directories = undetermined_in_root_folder
        else:
            no_sample_directories = False

        register_project_fastq_datafiles(
            run_id,
            proj_path,
            samplesheet,
            fq_dataset_url,
            uploader,
            fastqc_data=fqc_summary,
            no_sample_directories=no_sample_directories,
            fast_mode=options.fast)

    logger.info("Ingestion of run %s complete !", run_id)


@atexit.register
def _cleanup_tmp():
    global TMPDIRS
    for tmpdir in TMPDIRS:
        if exists(tmpdir):
            shutil.rmtree(tmpdir)
            logger.info("Removed temp directory: %s", tmpdir)


def run_in_console():
    MyTardisUploader.user_agent_name = os.path.basename(sys.argv[0])
    try:
        ingest_run()
    except Exception as e:
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


if __name__ == "__main__":
    run_in_console()
