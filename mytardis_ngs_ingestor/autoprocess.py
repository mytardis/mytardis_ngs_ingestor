import sys, os
from os import path
import subprocess
import jsondate as json
from datetime import datetime
import atexit
import toml
from toolz import dicttoolz

from argparse import ArgumentParser

import logging
from mytardis_uploader import setup_logging
logger = setup_logging()
from illumina import run_info
from illumina import fastqc
run_info.logger = logger
fastqc.logger = logger

# logger = logging.getLogger()

from simplekv.fs import FilesystemStore

from mytardis_ngs_ingestor.utils import config_helper
import mytardis_uploader
from illumina.run_info import get_run_id_from_path, get_sample_project_mapping
from illumina import fastqc
import illumina_uploader

# status enum
CREATED = 'created'    # record created, no action taken on task yet
RUNNING = 'running'
COMPLETE = 'complete'
ERROR = 'error'

current_task = None  # current task
taskdb = None


def first(l):
    """
    Return the first item from a list, or None of the list is empty.

    :param l: The list
    :type l: list
    :return: The first list item (or None)
    :rtype: object | None
    """
    return next(iter(l), None)


# def to_json(d):
#     """
#     Serialize a dictionary to JSON, correctly handling datetime.datetime
#     objects (to ISO 8601 dates, as strings).
#
#     :type d: dict
#     :rtype: str
#     """
#     if not isinstance(d, dict):
#         d = d.__dict__
#         # raise TypeError("Must be a dictionary")
#
#     date_handler = lambda obj: (
#         obj.isoformat(' ')
#         if isinstance(obj, datetime.date)
#         or isinstance(obj, datetime.datetime)
#         else None
#     )
#     return json.dumps(d, default=date_handler)


# TODO: alternative: use dataset (with sqlite) or jsonpickle or tinydb or lifter
class ProcessingTask:
    def __init__(self, run_id, task_name, status, timestamp=None, info=None):
        self.timestamp = timestamp or datetime.now()
        self.run_id = run_id
        self.task_name = task_name
        self.status = status  # running / complete / error
        self.info = info or {}

    # def as_dict(self):
    #     return {'timestamp': self.timestamp,
    #             'run_id': self.run_id,
    #             'task_name': self.task_name,
    #             'status': self.status,
    #             'info': self.info}
    #
    #     # return self.__dict__

    def is_complete(self):
        return self.status == COMPLETE

    def is_pending(self):
        return self.status != COMPLETE

    def is_failed(self):
        return self.status == ERROR

    def __str__(self):
        return json.dumps(self.__dict__) + '\n'


class TaskDb(FilesystemStore):
    def __init__(self, basedir, default_type, db_id=None, db_name='tasks'):
        self.dbpath = path.join(basedir, db_name)

        if not path.exists(self.dbpath):
            os.mkdir(path.join(basedir, db_name))

        super(TaskDb, self).__init__(self.dbpath)

        self.db_id = db_id
        self.default_type = default_type

    def _put(self, key, data):
        if not isinstance(data, dict):
            data = data.__dict__
        return super(TaskDb, self)._put(key, '%s\n' % json.dumps(data))

    def _get(self, key):
        value = super(TaskDb, self)._get(key)
        if value.strip():
            return json.loads(value)
        else:
            return None

    def get_as(self, key, *args, **kwargs):
        """
        Retrieve the object by it's key (as defined by TaskDb.default_type).
        Optionally, include a default return value if the key isn't found
        (otherwise raise KeyError). The keyword parameter 'object_type' can
        also be provided to override the type.

        eg.
          >> taskdb = TaskDb('/tmp', 'some_run_id', ProcessingTask)
          >> taskdb.get_as('nonexistanttask', 'a_default_return_value')
          'a_default_return_value'
          >>

        :param key:
        :type key:
        :param args:
        :type args:
        :param kwargs:
        :type kwargs:
        :return:
        :rtype:
        """
        object_type = kwargs.get('object_type', None)
        if object_type is None:
            object_type = self.default_type
        try:
            value = self.get(key)
            if value:
                return object_type(**value)
            else:
                return None
        except KeyError as e:
            if len(args) > 0:
                # default value
                return args[0]
            else:
                raise e

    def get_or_create(self, task_name):
        try:
            t = self.get_as(task_name)
        except KeyError as e:
            t = ProcessingTask(self.db_id, task_name, CREATED)
            self.put(task_name, t)
        return t

    def update(self, task):
        if isinstance(task, ProcessingTask):
            self.put(task.task_name, task)
        else:
            raise TypeError("task must be an instance of ProcessingTask")

    # def get_or_create(self, key, value):
    #     """
    #     Get value by key, with optional default value to return if missing.
    #
    #     :param value:
    #     :type value:
    #     :param key:
    #     :type key:
    #     :param args:
    #     :type args:
    #     :return:
    #     :rtype:
    #     """
    #     try:
    #         return self.get(key)
    #     except KeyError as e:
    #         self.put(key, value)
    #         return value

    def get_task_filepath(self, t):
        return path.join(self.dbpath, t.task_name)

    def is_complete(self, task_name):
        try:
            t = self.get_as(task_name)
            return t.status == COMPLETE
        except KeyError as e:
            return False

    def is_pending(self, task_name):
        return not self.is_complete(task_name)

    def exists(self, task_name):
        try:
            t = self.get_as(task_name)
            return t  # truthy
        except KeyError as e:
            return None  # falsey


def log_status(task, verbose=True):
    global taskdb

    log_fn = logger.info
    if task.status == ERROR:
        log_fn = logger.error

    if verbose:
        log_fn('run_id: %s task: %s status: %s',
               task.run_id.ljust(24),
               task.task_name.ljust(24),
               task.status)

        if task.status == ERROR and taskdb:
            log_fn('To retry this task, remove the file: %s',
                   taskdb.get_task_filepath(task))


# def run_task(task_name, fn, check_success_fn, verbose=True, *args, **kwargs):
#     current_task = taskdb.get_or_create(task_name)
#     if current_task.is_failed():
#         log_status(current_task, verbose)
#         return current_task
#
#     try:
#         current_task.status = RUNNING
#         taskdb.update(current_task)
#         if check_success_fn(*fn(*args, **kwargs)):
#             current_task.status = COMPLETE
#             taskdb.update(current_task)
#         else:
#             raise Exception('%s task failed.' % task_name)
#     except Exception as e:
#         logger.exception("%s task failed." % task_name)
#         current_task.status = ERROR
#         taskdb.update(current_task)
#
#     log_status(current_task, verbose)
#
#     if taskdb.is_pending(task_name):
#         return current_task


def do_create_checksum_manifest(run_dir, options):
    global current_task
    current_task = taskdb.get_or_create('create_checksum_manifest')

    if current_task.is_failed():
        log_status(current_task, options.verbose)
        return current_task

    cmd_out = None
    try:
        current_task.status = RUNNING
        # current_task.info['command'] = cmd
        taskdb.update(current_task)
        success, cmd_out = run_create_checksum_manifest(run_dir)
        current_task.status = COMPLETE
        current_task.info['output'] = cmd_out
        taskdb.update(current_task)
    except Exception as e:
        logger.exception("create_checksum_manifest task failed with an "
                         "exception.")
        current_task.status = ERROR
        current_task.info['output'] = cmd_out
        taskdb.update(current_task)

    log_status(current_task, options.verbose)

    return current_task


def do_rsync_to_archive_task(run_dir, options):
    task_name = 'rsync_to_archive'
    global current_task
    current_task = taskdb.get_or_create(task_name)

    # archive_basepath = '/srv/sonas/mhtp/market/illumina/'
    archive_basepath = '/data/illumina/archive_test'
    rsync_extra = '--exclude="{run_dir}/Thumbnail_Images/" ' \
                  '--include="*.xml" ' \
                  '--include="*.csv" ' \
                  '--include="*.log" ' \
                  '--exclude="{run_dir}/Data/Intensities/BaseCalls/L00*/" '.format(
        run_dir=run_dir)

    if current_task.is_failed():
        log_status(current_task, options.verbose)
        return current_task

    cmd_out = None
    try:
        current_task.status = RUNNING
        # current_task.info['command'] = cmd
        taskdb.update(current_task)
        success, cmd_out = run_rsync(run_dir, archive_basepath,
                                     sudo=True,
                                     chown='root:root',
                                     extra_args=rsync_extra)
        if success:
            current_task.status = COMPLETE
            current_task.info['output'] = cmd_out
            taskdb.update(current_task)
        else:
            raise Exception('%s task failed (non-zero exit code).' % task_name)
    except Exception as e:
        logger.exception("%s task failed." % task_name)
        current_task.status = ERROR
        current_task.info['output'] = cmd_out
        taskdb.update(current_task)

    log_status(current_task, options.verbose)

    return current_task


def try_autoprocessing(run_dir, options):
    global current_task
    global taskdb

    run_id = get_run_id_from_path(run_dir)
    taskdb = TaskDb(run_dir, ProcessingTask, run_id)

    ##
    # Check if the run should be ignored.
    # If so, silently skip it (unless we are being verbose)
    ##
    if taskdb.exists('ignore'):
        if options.verbose:
            logger.info("Skipping %s, set to ignore.", run_id)
        return taskdb.get_as('ignore')

    ##
    # Check if all processing is complete for the run.
    # If so, silently skip it (unless we are being verbose)
    ##
    current_task = taskdb.get_as('all_complete', None)
    if taskdb.is_complete('all_complete'):
        if options.verbose:
            log_status(current_task, options.verbose)
        return current_task

    if options.verbose:
        logger.info('Starting autoprocessing on: %s', run_dir)

    ##
    # Check if run has finished transferring from the sequencer
    ##
    current_task = taskdb.get_or_create('rta_complete')
    if current_task.is_pending():
        if run_is_complete(run_dir):
            current_task.status = COMPLETE
            taskdb.update(current_task)
        else:
            current_task.status = RUNNING
            taskdb.update(current_task)

    log_status(current_task, options.verbose)

    if taskdb.is_pending('rta_complete'):
        return current_task

    bcl2fastq_output_dir = path.join(run_dir, '%s.bcl2fastq' % run_id)
    # bcl2fastq_extra_args = ["--tiles s_1_1101"]
    bcl2fastq_extra_args = get_bcl2fastq_extra_args(run_dir)

    ##
    # Demultiplex with bcl2fastq
    ##
    current_task = taskdb.get_or_create('bcl2fastq')
    if current_task.is_failed():
        log_status(current_task, options.verbose)
        return current_task

    if current_task.is_pending():
        current_task.status = RUNNING
        taskdb.update(current_task)
        success, output = run_bcl2fastq(run_dir,
                                        output_directory=bcl2fastq_output_dir,
                                        stderr_file=path.join(run_dir, '%s.err' % run_id),
                                        extra_args=bcl2fastq_extra_args)
        if success:
            current_task.status = COMPLETE
            current_task.info = {'output': output}
            taskdb.update(current_task)
        else:
            current_task.status = ERROR
            current_task.info = {'output': output}
            taskdb.update(current_task)
            # log_status(current_task, options.verbose)
            # return current_task

    log_status(current_task, options.verbose)
    if current_task.is_pending():
        return current_task


    ##
    # Run FastQC on each Project in the run
    ##
    current_task = taskdb.get_or_create('fastqc')
    if current_task.is_failed():
        log_status(current_task, options.verbose)
        return current_task

    if current_task.is_pending():
        current_task.status = RUNNING
        taskdb.update(current_task)

        try:
            fastqs_per_project = get_sample_project_mapping(bcl2fastq_output_dir)
            ok = []
            failing_project = None
            for project, fastqs in fastqs_per_project.items():
                fastqs = [path.join(bcl2fastq_output_dir, fq) for fq in fastqs]
                proj_path = path.join(bcl2fastq_output_dir, project)
                result, output = fastqc.run_fastqc_on_project(fastqs,
                                                              proj_path,
                                                              clobber=True)
                ok.append(result)
                # fail early
                if result is None:
                    failing_project = project
                    break

            if all(ok):
                current_task.status = COMPLETE
                taskdb.update(current_task)
            else:
                current_task.status = ERROR
                current_task.info = {'project': failing_project}
                taskdb.update(current_task)
                log_status(current_task, options.verbose)
                return current_task
        except Exception as e:
            current_task.status = ERROR
            taskdb.update(current_task)
            log_status(current_task, options.verbose)
            logger.exception('FastQC task on %s raised an exception.' % run_dir)
            return current_task

    log_status(current_task, options.verbose)
    if current_task.is_pending():
        return current_task


    ##
    # Create checksum manifest file
    ##
    """
    current_task = do_create_checksum_manifest(run_dir, options)
    if current_task.is_failed() or current_task.is_pending():
        return current_task
    """

    ##
    # Custom scripts on run_dir / bcl2fastq_output_dir
    ##
    # TODO: eg, Copy to archival storage


    ##
    # Rsync to archival storage
    ##
    """
    current_task = do_rsync_to_archive_task(run_dir, options)
    if current_task.is_failed() or current_task.is_pending():
        return current_task
    """

    ##
    # Ingest using illumina_uploader
    ##
    current_task = taskdb.get_or_create('mytardis_upload')
    if not options.uploader:
        logger.exception("mytardis_upload task failed - config %s not found",
                         options.uploader_config)
        current_task.status = ERROR
        taskdb.update(current_task)

    if current_task.is_failed():
        log_status(current_task, options.verbose)
        return current_task

    try:
        current_task.status = RUNNING
        taskdb.update(current_task)
        illumina_uploader.logger = logger
        illumina_uploader.ingest_run(options.uploader, run_path=run_dir)
        current_task.status = COMPLETE
        taskdb.update(current_task)
    except Exception as e:
        logger.exception("mytardis_upload task failed with an exception.")
        current_task.status = ERROR
        taskdb.update(current_task)

    log_status(current_task, options.verbose)

    if taskdb.is_pending('mytardis_upload'):
        return current_task


    ##
    # Create an 'all_complete' task that allows us to silently skip
    # full processed runs.
    ##
    current_task = taskdb.get_or_create('all_complete')
    current_task.status = COMPLETE
    taskdb.update(current_task)
    log_status(current_task, options.verbose)

    logger.info('Autoprocessing completed for: %s', run_dir)

    return current_task


def is_illumina_run(run_dir):
    """
    Detects signature files in the run directory (eg RunInfo.xml) to detemine
    if it's likely to be an Illumina sequencing run or not.

    :param run_dir: The path to the run.
    :type run_dir: str
    :return: True if it looks like an Illumina run.
    :rtype: bool
    """

    # Ideas: Detecting RunInfo.xml is probably good enough, but we could
    #        look for other files, and analyse the directory name for the
    #        {date}_{instrument_id}_{run_number}_{flowcell_id} pattern too

    if path.exists(path.join(run_dir, 'RunInfo.xml')):
        return True
    return False


def run_is_complete(runfolder_dir, complete_file='RTAComplete.txt'):
    """
    Detect when a run is complete.

    :param runfolder_dir:
    :type runfolder_dir:
    :param complete_file:
    :type complete_file:
    :return:
    :rtype:
    """
    return path.exists(path.join(runfolder_dir, complete_file))


def run_bcl2fastq(runfolder_dir,
                  output_directory=None,
                  bcl2fastqc_bin=None,
                  stderr_file='',
                  nice=False,
                  version='2.19',
                  extra_args=None):
    """
    Run bcl2fastq with commandline options.

    :param runfolder_dir:
    :type runfolder_dir:
    :param output_directory:
    :type output_directory:
    :param bcl2fastqc_bin:
    :type bcl2fastqc_bin:
    :param nice:
    :type nice:
    :param extra_args:
    :type extra_args:
    :return: A tuple of (output_directory, stdout+stderr),
             either of which may be None
    :rtype: tuple(str | None, str | None)
    """

    if not runfolder_dir or \
            not isinstance(runfolder_dir, basestring) or \
            not path.exists(runfolder_dir):
        raise ValueError("runfolder_dir must be a valid path string")

    if not bcl2fastqc_bin:
        # We will assume it's on path with the default name
        bcl2fastqc_bin = 'bcl2fastq'

    if nice:
        nice = 'nice '
    else:
        nice = ''

    options = ['--runfolder-dir %s' % runfolder_dir]

    if output_directory:
        options.append('--output-dir %s' % output_directory)

    if extra_args:
        options.extend(extra_args)

    if stderr_file:
        stderr_file = ' 2>&1 | tee -a %s' % stderr_file

    cmd = '{nice}{bcl2fastq} {options} {stderr_file}'\
        .format(nice=nice,
                bcl2fastq=bcl2fastqc_bin,
                stderr_file=stderr_file,
                options=' '.join(options))

    # eg, using a Docker container prepared like:
    # https://gist.github.com/pansapiens/0e9b36cc1b11ce3c6e49dc81d09e30bf
    # cmd = '{nice} docker run -it ' \
    #       '-v {output_directory}:/output ' \
    #       '-v {runfolder_dir}:/run bcl2fastq:{version} ' \
    #       '/usr/local/bin/bcl2fastq -o /output -R /run'.format(
    #     nice=nice,
    #     version=version,
    #     output_directory=output_directory,
    #     runfolder_dir=runfolder_dir,
    # )

    logger.info('Running bcl2fastq on: %s', runfolder_dir)

    cmd_out = None
    try:
        # bcl2fastq writes everything useful to stderr
        cmd_out = subprocess.check_output(cmd,
                                          shell=True,
                                          stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        logger.error('bcl2fastq failed: stdout/stderr: %s', cmd_out)
        return None, cmd_out

    success = False
    # A successful outcome has something like this as the last line in stderr
    # "2017-04-02 15:56:32 [1bbb880] Processing completed with 0 errors and 0 warnings."
    if cmd_out and 'Processing completed with 0 errors' in cmd_out.splitlines()[-1]:
        success = True
    else:
        logger.error('bcl2fastq failed stdout/stderr: %s', cmd_out)

    if not success:
        return None, cmd_out

    return output_directory, cmd_out


def get_bcl2fastq_extra_args(run_dir, filename='bcl2fastq_extra_args.txt'):
    extra_args_path = path.join(run_dir, filename)
    lines = []
    if path.isfile(extra_args_path):
        with open(extra_args_path, 'r') as f:
            lines = f.readlines()
    lines = [l.replace('\n', ' ').replace('\r', ' ') for l in lines]
    return lines


def run_create_checksum_manifest(base_dir,
                                 manifest_filename='manifest-%s.txt',
                                 hash_type='md5'):


    manifest_filename %= hash_type
    manifest_filepath = path.join(base_dir, manifest_filename)
    cmd = 'find {base_dir} ' \
          '-type f ' \
          '-exec {hash_type}sum "{{}}" + >{manifest_filepath}'.format(
        base_dir=base_dir,
        hash_type=hash_type,
        manifest_filepath=manifest_filepath)

    cmd_out = None
    try:
        cmd_out = subprocess.check_output(cmd,
                                          shell=True,
                                          stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        logger.error('create checksum failed: stdout/stderr: %s', cmd_out)
        return None, cmd_out

    return manifest_filepath, cmd_out


def run_rsync(src, dst, sudo=False, chown=None, checksum=True, extra_args=''):

    if checksum:
        extra_args += ' --checksum '
    if chown:
        extra_args += ' --chown=%s ' % chown
    if sudo:
        sudo = 'sudo'
    else:
        sudo = ''

    cmd = '{sudo} rsync -av ' \
          ' {extra_args} ' \
          '{src} {dst}/'.format(sudo=sudo,
                                 src=src,
                                 dst=dst,
                                 extra_args=extra_args)

    cmd_out = None
    try:
        cmd_out = subprocess.check_output(cmd,
                                          shell=True,
                                          stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        logger.error('rsync failed: stdout/stderr: %s', cmd_out)
        return None, cmd_out

    return dst, cmd_out


def process_all_runs(run_storage_base, options):
    """

    Returns True if there were no error, False otherwise.

    :param run_storage_base:
    :type run_storage_base:
    :param options:
    :type options:
    :return:
    :rtype:
    """
    global logging
    global current_task

    errored_tasks = []
    for name in os.listdir(run_storage_base):
        run_dir = path.join(run_storage_base, name)
        run_id = get_run_id_from_path(run_dir)
        if path.isdir(run_dir) and is_illumina_run(run_dir):
            try:
                emitted_task = try_autoprocessing(run_dir, options)
                if emitted_task and emitted_task.status == ERROR:
                    errored_tasks.append(emitted_task)
            except Exception as e:
                # Dummy catchall task to signal exceptional failure
                errored_tasks.append(ProcessingTask(run_id,
                                                    'try_autoprocessing',
                                                    ERROR))
            current_task = None

            # Stop on error in any task, don't continue with the other runs
            # if emitted_task and emitted_task.status != COMPLETE:
            #     break

    if options.verbose:
        if errored_tasks:
            errorlist = ', '.join(['%s:%s' % (t.task_name, t.run_id)
                                 for t in errored_tasks])
            logger.error("Processing runs in %s completed with failures: %s",
                          run_storage_base, errorlist)
        else:
            logger.info("Successfully completed processing of runs in: %s",
                         run_storage_base)

    return not errored_tasks


def watch_runs(run_storage_base, options, tight_loop_timer=1):
    """
    Typically this script will be run via cron, however if for some
    reason you'd rather run it continuously (under supervisor) this will do it.

    :param tight_loop_timer:
    :type tight_loop_timer:
    :return:
    :rtype:
    """
    import time
    import crython

    # @crython.job(minute=range(0, 60, frequency), second=0)
    # @crython.job(second=[5])  # once every minute
    # @crython.job(second=range(0, 60, 5))  # once every 5 seconds
    @crython.job(minute=range(0, 60, 15))  # once every 15 mins
    def _process_all_cron():
        logger.info('Running scheduled autoprocessing on: %s', run_storage_base)
        process_all_runs(run_storage_base, options)

    crython.start()

    try:
        while True:
            time.sleep(tight_loop_timer)
    except KeyboardInterrupt:
        logging.info("Stopped watching: %s" % run_storage_base)
        pass


def get_commandline_args(parser=None):
    """
    Parses commandline options (extending an existing parser),
    returns argparse options object.

    :param parser: The parser to extend.
    :type parser: argparse.ArgumentParser
    :return: The extended parser.
    :rtype: argparse.ArgumentParser
    """
    if not parser:
        parser = ArgumentParser()

    parser.add_argument("--quiet",
                        action="store_false",
                        dest="verbose",
                        default=True,
                        help="Less verbose logging. When set, a subsequent "
                             "successful walk over a set of processed runs "
                             "where not additional processing occurred will"
                             "be silent.")

    # parser.add_argument('--config',
    #                     dest="config_file",
    #                     type=str,
    #                     metavar="AUTOPROCESSNG_CONFIG")

    parser.add_argument("--uploader-config",
                        dest="uploader_config",
                        type=str,
                        default="uploader_config.toml",
                        help="The path to the uploader config file "
                             "eg uploader_config.toml",
                        metavar="UPLOADER_CONFIG")

    parser.add_argument("--runs",
                        dest="run_storage_base",
                        type=str,
                        help="The absolute PATH to a directory containing "
                             "multiple runs to be processed (eg "
                             "/data/illumina)",
                        metavar="RUNS_STORAGE_BASE")

    parser.add_argument("--single-run",
                        dest="run_path",
                        type=str,
                        help="The absolute PATH to a single run to be "
                             "processed (eg "
                             "/data/illumina/170404_SNL177_0169_AHHGVYBCXY)")

    parser.add_argument("--watch",
                        action="store_true",
                        dest="watch",
                        default=False,
                        help="An alternative to running under cron - remain "
                             "running and watch for new runs. "
                             "Stop with Ctrl-C.")

    # parser.add_argument("-r", "--dry",
    #                     action="store_true",
    #                     dest="dry_run",
    #                     default=False,
    #                     help="Dry run (don't actually process, just show "
    #                          "what the next task would be)")

    # options = parser.parse_args()
    # return parser, options
    return parser


@atexit.register
def _set_current_task_status_on_exit():
    """
    Attempt to set current task status to 'error' if we have a premature
    exit.
    """
    global current_task
    global taskdb

    if current_task is not None \
            and taskdb is not None \
            and current_task.status != COMPLETE:
        current_task.status = ERROR
        taskdb.put(current_task.task_name, current_task)


def run_in_console():
    global logger

    parser = get_commandline_args()
    options = parser.parse_args()
    # options, argv = parser.parse_known_args()

    options.uploader = None
    if options.uploader_config and os.path.isfile(options.uploader_config):
        uploader_parser = ArgumentParser()
        uploader_parser = mytardis_uploader.setup_commandline_args(
            uploader_parser)
        uploader_parser = illumina_uploader.setup_commandline_args(
            uploader_parser)

        uploader_default_options, argv = uploader_parser.parse_known_args()
        # uploader_default_options = uploader_parser.parse_args()

        try:
            uploader_options = config_helper.get_config_toml(
                options.uploader_config)
        except IOError as e:
            parser.error("Cannot read config file: %s" %
                         options.uploader_config)
            sys.exit(1)

        # Any option that is unset in the commandline args will
        # receive a value from a config file option
        uploader_parser.set_defaults(**uploader_options)
        # TODO: By using parse_known_args rather than parse_args, we might
        #       inadvertedly allow typos in config file keys to be ignored ?
        options.uploader, argv = uploader_parser.parse_known_args()
        # options.uploader = uploader_parser.parse_args()

    # options.uploader = uploader_options

    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit()
    try:
        if options.run_storage_base and options.run_path:
            logger.error("Please use only --runs or --single-run, not both.")

        if options.watch and not options.run_storage_base:
            logger.error("You must specify the --runs option if using --watch.")

        if options.run_storage_base:
            if options.watch:
                watch_runs(options.run_storage_base, options)
                sys.exit()
            else:
                logger.info('Running single autoprocess pass on: %s',
                            options.run_storage_base)
                exit_code = 0
                ok = process_all_runs(options.run_storage_base, options)
                if not ok:
                    exit_code = 1
                sys.exit(exit_code)

        if options.run_path:
            try_autoprocessing(options.run_path, options)
            sys.exit()
    except KeyboardInterrupt as e:
        _set_current_task_status_on_exit()


if __name__ == '__main__':
    run_in_console()
