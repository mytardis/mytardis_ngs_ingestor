import os
from os import path
import subprocess
from datetime import datetime
import jsondate as json

from simplekv.fs import FilesystemStore

import logging

from mytardis_ngs_ingestor.illumina import run_info

from mytardis_ngs_ingestor.illumina.run_info import get_run_id_from_path, \
    get_sample_project_mapping
from mytardis_ngs_ingestor.illumina import fastqc
from mytardis_ngs_ingestor import illumina_uploader

# status enum
CREATED = 'created'  # record created, no action taken on task yet
RUNNING = 'running'
COMPLETE = 'complete'
ERROR = 'error'


class Current(object):
    _singleton = None

    def __init__(self):
        self.task = None

    def __new__(cls, *args, **kwargs):
        if not cls._singleton:
            cls._singleton = super(Current, cls).__new__(cls, *args, **kwargs)
        return cls._singleton


class ProcessingTask(object):
    def __init__(self, run_id, task_name, status,
                 db=None,
                 timestamp=None,
                 last_failure_notify_time=None,
                 info=None):
        self.timestamp = timestamp or datetime.now()
        self.last_failure_notify_time = last_failure_notify_time
        self.run_id = run_id
        self.task_name = task_name
        self.status = status  # running / complete / error
        self.info = info or {}
        self._db = db  # underscore = not serialized

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

    def is_running(self):
        return self.status == RUNNING

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
        import jsondate as json

        if not isinstance(data, dict):
            data = dict(data.__dict__)
            # any key starting with underscore isn't serialized
            for k in data.keys():
                if k.startswith('_'):
                    del data[k]

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
            t = ProcessingTask(self.db_id, task_name, CREATED, db=self)
            self.put(task_name, t)
        return t

    def update(self, task):
        if isinstance(task, ProcessingTask):
            self.put(task.task_name, task)
        else:
            raise TypeError("task must be an instance of ProcessingTask")

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
    log_fn = logging.info
    if task.status == ERROR:
        log_fn = logging.error

    if verbose:
        log_fn('run_id: %s task: %s status: %s',
               task.run_id.ljust(24),
               task.task_name.ljust(24),
               task.status)

        if task.status == ERROR:
            task_path = task.task_name
            if task._db:
                task_path = task._db.get_task_filepath(task)
            log_fn('To retry this task, remove the file: %s',
                   task_path)


def log_retry(task, verbose=True):
    if verbose:
        logging.warning('Retrying - run_id: %s task: %s',
                        task.run_id.ljust(24),
                        task.task_name.ljust(24))


def run_bcl2fastq(runfolder_dir,
                  output_directory=None,
                  bcl2fastq_bin=None,
                  stderr_file='',
                  nice=False,
                  docker_image=None,
                  extra_args=None):
    """
    Run bcl2fastq with commandline options.
    If docker_image is specified, run inside a Docker container.

    :param stderr_file:
    :type stderr_file:
    :param docker_image: A Docker image ang tag (eg
                         "genomicpariscentre/bcl2fastq2:2.19.0.316"). If
                         specified, bcl2fastq is run in a container.
    :type docker_image: str
    :param runfolder_dir:
    :type runfolder_dir:
    :param output_directory:
    :type output_directory:
    :param bcl2fastq_bin:
    :type bcl2fastq_bin:
    :param nice:
    :type nice:
    :param extra_args:
    :type extra_args:
    :return: A tuple of (success, stdout+stderr)
    :rtype: tuple(bool, str | None)
    """

    if not runfolder_dir or \
            not isinstance(runfolder_dir, basestring) or \
            not path.exists(runfolder_dir):
        raise ValueError("runfolder_dir must be a valid path string")

    if not bcl2fastq_bin:
        # We will assume it's on path with the default name
        bcl2fastq_bin = 'bcl2fastq'

    if nice:
        nice = 'nice '
    else:
        nice = ''

    options = []

    if extra_args:
        options.extend(extra_args)

    if stderr_file:
        stderr_file = ' 2>&1 | tee -a %s' % stderr_file

    if docker_image is not None:
        if not output_directory:
            output_directory = path.join(
                runfolder_dir,
                "%s.bcl2fastq" % path.basename(runfolder_dir))

        # eg, using a Docker container prepared like:
        # https://gist.github.com/pansapiens/0e9b36cc1b11ce3c6e49dc81d09e30bf
        cmd = 'docker run -it ' \
              '-v {output_directory}:/output ' \
              '-v {runfolder_dir}:/run {docker_image} ' \
              '{nice} {bcl2fastq} ' \
              '--output-dir /output ' \
              '--runfolder-dir /run ' \
              '{options} {stderr_file}'.format(
                nice=nice,
                docker_image=docker_image,
                bcl2fastq=bcl2fastq_bin,
                output_directory=output_directory,
                runfolder_dir=runfolder_dir,
                options=' '.join(options),
                stderr_file=stderr_file)

        # If the current user/group (by UID/GID) exists inside the container
        # you can use this
        # '--user `id -n -u`:`id -n -g` '

    else:
        options.append('--runfolder-dir %s' % runfolder_dir)

        if output_directory:
            options.append('--output-dir %s' % output_directory)

        cmd = '{nice}{bcl2fastq} {options} {stderr_file}'.format(
            nice=nice,
            bcl2fastq=bcl2fastq_bin,
            stderr_file=stderr_file,
            options=' '.join(options))

    logging.info('Running bcl2fastq on: %s', runfolder_dir)
    logging.info('Command: %s', cmd)

    cmd_out = None
    success = False
    try:
        # bcl2fastq writes everything useful to stderr
        cmd_out = subprocess.check_output(cmd,
                                          shell=True,
                                          stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        logging.error('bcl2fastq failed: stdout/stderr: %s', cmd_out)
        return success, cmd_out

    # We don't trust error codes in this case, we look at the stderr messages.
    # A successful outcome has something like this as the last line in stderr
    # "2017-04-02 15:56:32 [1bbb880] Processing completed with 0 errors and 0 warnings."
    if cmd_out and 'Processing completed with 0 errors' in cmd_out:
        success = True
    else:
        logging.error('bcl2fastq failed stdout/stderr: %s', cmd_out)

    return success, cmd_out


def get_legcay_bcl2fastq_extra_args(run_dir,
                                    filename='bcl2fastq_extra_args.txt'):
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
    # TODO: Make an equivalent task that uses hashdeep instead
    #       (you can sudo apt install hashdeep in 16.04)
    #       Or optionally write hashdeep format here
    #       https://gist.github.com/techtonik/5175896

    manifest_filename %= hash_type
    manifest_filepath = path.join(base_dir, manifest_filename)
    cmd = 'find {base_dir} ' \
          '-type f ' \
          '-exec {hash_type}sum "{{}}" + >{manifest_filepath}'.format(
            base_dir=base_dir,
            hash_type=hash_type,
            manifest_filepath=manifest_filepath)

    logging.info('Command: %s', cmd)

    cmd_out = None
    try:
        cmd_out = subprocess.check_output(cmd,
                                          shell=True,
                                          stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        logging.error('create checksum failed: stdout/stderr: %s', cmd_out)
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

    logging.info("Command: %s", cmd)

    cmd_out = None
    try:
        cmd_out = subprocess.check_output(cmd,
                                          shell=True,
                                          stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        logging.error('rsync failed: stdout/stderr: %s', cmd_out)
        return None, cmd_out

    return dst, cmd_out


def do_rta_complete(taskdb, current, run_dir, options):
    current.task = taskdb.get_or_create('rta_complete')

    if options.retry and current.task.is_failed():
        current.task.status = CREATED
        log_retry(current.task, options.verbose)

    # status == RUNNING for this task means a run directory exists but
    # RTAComplete.txt doesn't.
    # So unlike most tasks we still need to allow the is_pending
    # check below to proceed even when the status == RUNNING, otherwise
    # status can never become COMPLETE;
    if current.task.is_failed():  # or current.task.is_running():
        log_status(current.task, options.verbose)
        return current.task

    if current.task.is_pending():
        if run_info.run_is_complete(run_dir):
            current.task.status = COMPLETE
            taskdb.update(current.task)
        else:
            current.task.status = RUNNING
            taskdb.update(current.task)

    log_status(current.task, options.verbose)
    return current.task


def do_bcl2fastq(taskdb, current, run_dir, options):
    current.task = taskdb.get_or_create('bcl2fastq')

    if options.retry and current.task.is_failed():
        current.task.status = CREATED
        log_retry(current.task, options.verbose)

    if current.task.is_failed() or current.task.is_running():
        log_status(current.task, options.verbose)
        return current.task

    if current.task.is_pending():
        current.task.status = RUNNING
        taskdb.update(current.task)

        opts = options.config.get('bcl2fastq', {})

        # output_dir = path.join(run_dir, '%s.bcl2fastq' % options.run_id)
        output_dir = None
        output_dir_fmt_str = opts.get('output-dir', None)
        if output_dir_fmt_str:
            output_dir = illumina_uploader.get_bcl2fastq_output_dir(
                output_dir_fmt_str,
                options.run_id,
                run_dir)

        config_args = opts.get('args', {})
        extra_args = []
        for arg, value in config_args.items():
            if isinstance(value, bool):
                value = ''
            extra_args.append("--%s %s" % (arg, value))

        success, output = run_bcl2fastq(
            run_dir,
            bcl2fastq_bin=opts.get('binary_path', None),
            docker_image=opts.get('docker_image', None),
            output_directory=output_dir,
            stderr_file=path.join(run_dir, '%s.err' % options.run_id),
            extra_args=extra_args)
        if success:
            current.task.status = COMPLETE
            current.task.info = {'output': output}
            taskdb.update(current.task)
        else:
            current.task.status = ERROR
            current.task.info = {'output': output}
            taskdb.update(current.task)
            # log_status(current.task, options.verbose)
            # return current.task

    log_status(current.task, options.verbose)
    return current.task


def do_fastqc(taskdb, current, run_dir, options):
    current.task = taskdb.get_or_create('fastqc')

    if options.retry and current.task.is_failed():
        current.task.status = CREATED
        log_retry(current.task, options.verbose)

    if current.task.is_failed() or current.task.is_running():
        log_status(current.task, options.verbose)
        return current.task

    if current.task.is_pending():
        current.task.status = RUNNING
        taskdb.update(current.task)

        try:
            outdir = '%s/Data/Intensities/BaseCalls/' % run_dir
            bcl2fastq_opts = options.config.get('bcl2fastq', {})
            outdir_fmt_str = bcl2fastq_opts.get('output-dir', None)
            if outdir_fmt_str:
                outdir = illumina_uploader.get_bcl2fastq_output_dir(
                    outdir_fmt_str,
                    options.run_id,
                    run_dir)

            fastqc_opts = options.config.get('fastqc', {})
            fastqc_bin = fastqc_opts.get('binary_path', None)
            fastqs_per_project = get_sample_project_mapping(outdir)
            ok = []
            failing_project = None
            for project, fastqs in fastqs_per_project.items():
                fastqs = [path.join(outdir, fq)
                          for fq in fastqs]
                proj_path = path.join(outdir, project)
                result, output = fastqc.run_fastqc_on_project(
                    fastqs,
                    proj_path,
                    fastqc_bin=fastqc_bin,
                    clobber=True)

                if 'Failed to process' in output:
                    result = None
                ok.append(result)
                # fail early
                if result is None:
                    failing_project = project
                    break

            if all(ok):
                current.task.status = COMPLETE
                taskdb.update(current.task)
            else:
                current.task.status = ERROR
                current.task.info = {'project': failing_project}
                taskdb.update(current.task)

        except Exception as e:
            current.task.status = ERROR
            taskdb.update(current.task)
            log_status(current.task, options.verbose)
            logging.exception(
                'FastQC task on %s raised an exception.' % run_dir)
            return current.task

    log_status(current.task, options.verbose)
    return current.task


def do_create_checksum_manifest(taskdb, current, run_dir, options):
    current.task = taskdb.get_or_create('create_checksum_manifest')

    if options.retry and current.task.is_failed():
        current.task.status = CREATED
        log_retry(current.task, options.verbose)

    if current.task.is_failed() or current.task.is_running():
        log_status(current.task, options.verbose)
        return current.task

    if current.task.is_pending():
        cmd_out = None
        try:
            current.task.status = RUNNING
            # current.task.info['command'] = cmd
            taskdb.update(current.task)
            success, cmd_out = run_create_checksum_manifest(run_dir)
            current.task.status = COMPLETE
            current.task.info['output'] = cmd_out
            taskdb.update(current.task)
        except Exception as e:
            logging.exception("create_checksum_manifest task failed with an "
                              "exception.")
            current.task.status = ERROR
            current.task.info['output'] = cmd_out
            taskdb.update(current.task)

    log_status(current.task, options.verbose)
    return current.task


def do_rsync_to_archive(taskdb, current, run_dir, options):
    task_name = 'rsync_to_archive'
    current.task = taskdb.get_or_create(task_name)

    if options.retry and current.task.is_failed():
        current.task.status = CREATED
        log_retry(current.task, options.verbose)

    if current.task.is_failed() or current.task.is_running():
        log_status(current.task, options.verbose)
        return current.task

    opts = options.config.get(task_name, {})
    target_basepath = opts.get('target_basepath', None)

    if not target_basepath:
        logging.exception("%s task failed - target_basepath not specified" %
                          task_name)
        current.task.status = ERROR
        taskdb.update(current.task)

        log_status(current.task, options.verbose)
        return current.task

    rsync_extra = opts.get('args', [])
    extra_args = ' '.join(rsync_extra)
    extra_args = extra_args.format(run_dir=run_dir)
    sudo = opts.get('sudo', False)
    chown = opts.get('chown', None)

    cmd_out = None
    try:
        current.task.status = RUNNING
        # current.task.info['command'] = cmd
        taskdb.update(current.task)
        success, cmd_out = run_rsync(run_dir, target_basepath,
                                     sudo=sudo,
                                     chown=chown,
                                     extra_args=extra_args)
        if success:
            current.task.status = COMPLETE
            current.task.info['output'] = cmd_out
            taskdb.update(current.task)
        else:
            raise Exception('%s task failed (non-zero exit code).' % task_name)
    except Exception as e:
        logging.exception("%s task failed." % task_name)
        current.task.status = ERROR
        current.task.info['output'] = cmd_out
        taskdb.update(current.task)

    log_status(current.task, options.verbose)
    return current.task


def do_mytardis_upload(taskdb, current, run_dir, options):
    current.task = taskdb.get_or_create('mytardis_upload')

    if options.retry and current.task.is_failed():
        current.task.status = CREATED
        log_retry(current.task, options.verbose)

    if current.task.is_failed() or current.task.is_running():
        log_status(current.task, options.verbose)
        return current.task

    if not options.config.get('mytardis_uploader', None):
        logging.exception("mytardis_upload task failed - config file %s not "
                          "found",
                          options.uploader_config)
        current.task.status = ERROR
        taskdb.update(current.task)
        return current.task

    try:
        current.task.status = RUNNING
        taskdb.update(current.task)
        options.config['mytardis_uploader']['path'] = run_dir
        illumina_uploader.ingest_run(options.config.mytardis_uploader,
                                     run_path=run_dir)
        current.task.status = COMPLETE
        taskdb.update(current.task)
    except Exception as e:
        logging.exception("mytardis_upload task failed with an exception.")
        current.task.status = ERROR
        taskdb.update(current.task)

    log_status(current.task, options.verbose)
    return current.task

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


# def run_task(task_name, fn, check_success_fn, verbose=True, *args, **kwargs):
#     current.task = taskdb.get_or_create(task_name)
#     if current.task.is_failed():
#         log_status(current.task, verbose)
#         return current.task
#
#     try:
#         current.task.status = RUNNING
#         taskdb.update(current.task)
#         if check_success_fn(*fn(*args, **kwargs)):
#             current.task.status = COMPLETE
#             taskdb.update(current.task)
#         else:
#             raise Exception('%s task failed.' % task_name)
#     except Exception as e:
#         logging.exception("%s task failed." % task_name)
#         current.task.status = ERROR
#         taskdb.update(current.task)
#
#     log_status(current.task, verbose)
#
#     if taskdb.is_pending(task_name):
#         return current.task
