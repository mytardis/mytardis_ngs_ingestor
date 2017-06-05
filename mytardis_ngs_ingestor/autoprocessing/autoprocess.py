import sys, os
import logging
import types
from os import path
from datetime import datetime, timedelta
import atexit
from attrdict import AttrDict
import collections

from argparse import ArgumentParser

from mytardis_ngs_ingestor.mytardis_uploader import setup_logging
from mytardis_ngs_ingestor import illumina_uploader
from mytardis_ngs_ingestor.illumina import run_info
from mytardis_ngs_ingestor.illumina import fastqc
from mytardis_ngs_ingestor.utils import config_helper
from mytardis_ngs_ingestor.illumina.run_info import (get_run_id_from_path,
                                                     get_sample_project_mapping)
import tasks
from tasks import TaskDb, ProcessingTask

current = tasks.Current()  # current task 'pseudo-singleton'
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


def get_task_functions(tasknames):
    """
    Given a list of task names, return a list of corresponding functions.
    Raises ValueError if an approriately named function can't be found.

    :param tasknames: A list of task names.
    :type tasknames: list(str)
    :return: A list of functions.
    :rtype: list(type.FunctionType)
    """
    task_fns = []
    for task_name in tasknames:
        # task_func = globals().get('do_%s' % task_name, None)
        task_func = getattr(tasks, 'do_%s' % task_name, None)
        if not task_func or not type(task_func) is types.FunctionType:
            logging.error("Task name '%s' is invalid.", task_name)
            raise ValueError("Task name '%s' is invalid.", task_name)
        task_fns.append(task_func)

    return task_fns


def try_autoprocessing(run_dir, options):
    global current
    global taskdb

    options = parse_run_specific_config(run_dir, options)

    run_id = get_run_id_from_path(run_dir)
    options.run_id = run_id

    taskdb = TaskDb(run_dir, ProcessingTask, run_id)

    def can_readwrite(p):
        return os.access(p, os.R_OK) and os.access(p, os.W_OK)

    if options.skip_bad_permissions:
        if not can_readwrite(run_dir) or not can_readwrite(taskdb.dbpath):
            logging.info("Skipping %s, no read/write permissions.", run_dir)
            return None

    ##
    # Check if the run should be ignored.
    # If so, silently skip it (unless we are being verbose)
    ##
    # if taskdb.exists('ignore'):
    if path.exists(path.join(taskdb.dbpath, 'ignore')):
        if options.verbose:
            logging.info("Skipping %s, set to ignore.", run_id)
        return taskdb.get_as('ignore')

    ##
    # Check if all processing is complete for the run.
    # If so, silently skip it (unless we are being verbose)
    ##
    current.task = taskdb.get_as('all_complete', None)
    if taskdb.is_complete('all_complete'):
        if options.verbose:
            tasks.log_status(current.task, options.verbose)
        return current.task

    if options.verbose:
        logging.info('Starting autoprocessing on: %s', run_dir)

    # Only execute tasks in options.config.tasks
    # Pre-check that all task names map to valid functions
    task_fns = []
    try:
        task_fns = get_task_functions(options.config.get('tasks', []))
    except ValueError as e:
        raise e

    if not task_fns:
        logging.warning("No tasks specified in config.")

    # Now actually run the tasks
    for task_fn in task_fns:
        current.task = task_fn(taskdb, current, run_dir, options)
        if current.task.is_failed() or current.task.is_pending():
            return current.task

    ##
    # Create an 'all_complete' task that allows us to silently skip
    # full processed runs.
    ##
    current.task = taskdb.get_or_create('all_complete')
    current.task.status = tasks.COMPLETE
    taskdb.update(current.task)
    tasks.log_status(current.task, options.verbose)

    done_msg = 'Autoprocessing completed for: %s' % run_dir
    logging.info(done_msg)
    logging.getLogger('autoprocess_notify').info(done_msg)

    return current.task


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
    global current

    errored_tasks = []
    running_tasks = []
    for name in os.listdir(run_storage_base):
        run_dir = path.join(run_storage_base, name)
        run_id = get_run_id_from_path(run_dir)
        if path.isdir(run_dir) and run_info.is_illumina_run(run_dir):
            try:
                emitted_task = try_autoprocessing(run_dir, options)
                if emitted_task and emitted_task.status == tasks.ERROR:
                    errored_tasks.append(emitted_task)
                if emitted_task and emitted_task.status == tasks.RUNNING:
                    running_tasks.append(emitted_task)
            except Exception as e:
                # Dummy catchall task to signal exceptional failure
                errored_tasks.append(ProcessingTask(run_id,
                                                    'try_autoprocessing',
                                                    tasks.ERROR))
                if options.verbose:
                    logging.exception(e)

            current.task = None

            # Stop on error in any task, don't continue with the other runs
            # if emitted_task and emitted_task.status != COMPLETE:
            #     break

    errorlist = ', '.join(['%s:%s' % (t.task_name, t.run_id)
                           for t in errored_tasks])
    runninglist = ', '.join(['%s:%s' % (t.task_name, t.run_id)
                             for t in running_tasks])
    if options.verbose:
        if running_tasks:
            logging.info("%s task(s) are currently running (%s): %s",
                         len(running_tasks), run_storage_base, runninglist)
        else:
            logging.info("Successfully completed processing of runs in: %s",
                         run_storage_base)

    if errored_tasks:
        if options.verbose:
            logging.error("Processing runs in %s completed with failures: %s",
                          run_storage_base, errorlist)
            logging.getLogger('autoprocess_notify').error(
                "Processing runs in %s completed with failures: %s",
                run_storage_base,
                errorlist)
        else:
            # Throttle notification log if in --quiet mode
            notify_every = timedelta(minutes=options.notify_frequency)
            for t in errored_tasks:
                if t.last_failure_notify_time is None or \
                        (t.last_failure_notify_time +
                            notify_every < datetime.now()):
                    logging.getLogger('autoprocess_notify').error(
                        "Processing runs in %s completed with failures: %s",
                        run_storage_base,
                        errorlist)
                    t.last_failure_notify_time = datetime.now()
                    taskdb = TaskDb(path.join(run_storage_base, t.run_id),
                                    ProcessingTask, t.run_id)
                    taskdb.update(t)
                    # t._db.update(t)
                    break

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
        logging.info('Running scheduled autoprocessing on: %s', run_storage_base)
        process_all_runs(run_storage_base, options)

    crython.start()

    try:
        while True:
            time.sleep(tight_loop_timer)
    except KeyboardInterrupt:
        logging.info("Stopped watching: %s" % run_storage_base)
        pass


def setup_commandline_args(parser=None):
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

    parser = _add_uploader_config_argparser(parser=parser)

    parser.add_argument("--quiet",
                        action="store_false",
                        dest="verbose",
                        # default=True,  # don't use
                        help="Less verbose logging. When set, a subsequent "
                             "successful walk over a set of processed runs "
                             "where not additional processing occurred will"
                             "be silent.")

    parser.add_argument("--retry",
                        action="store_true",
                        dest="retry",
                        help="Removes any failed tasks from previous "
                             "invocations and allows them to be retried.")

    parser.add_argument('--config',
                        dest="config_file",
                        type=str,
                        default='autoprocess_config.toml',
                        help="The global config file to use for autoprocessing "
                             "settings. A config file "
                             "'autoprocessing_config.toml' in individual run "
                             "directories overrides settings in this file."
                             "Commandline options override all config file"
                             "settings.",
                        metavar="AUTOPROCESSNG_CONFIG")

    parser.add_argument("--logging-config",
                        dest="logging_config",
                        type=str,
                        # default="logging_config.toml",
                        help="The path to the logging config file "
                             "eg logging_config.toml",
                        metavar="LOGGING_CONFIG")

    # TODO: It might be better to make these subparser modes like:
    #       autoprocess process --run-dir /data/runs
    #       autoprocess process --single-run /data/runs/blabla
    #       autoprocess watch --run-dir /data/runs
    #       # Wait for a single run to become complete, process then exit
    #       autoprocess watch --single-run /data/runs/blabla
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
                        # default=False, # don't use
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


def _add_uploader_config_argparser(parser=None):
    """
    We use this to add the --uploader-config argument to an ArgParser instance.
    It exists as a way to create an ArguementParser that ONLY recognises this
    option, to be passed to illumina_uploader.get_config_options.

    Feel like an ugly hack. Sorry :/

    :param parser: The parser to extend.
    :type parser: argparse.ArgumentParser
    :return: The extended parser.
    :rtype: argparse.ArgumentParser
    """
    if not parser:
        parser = ArgumentParser()

    parser.add_argument("--uploader-config",
                        dest="uploader_config",
                        type=str,
                        # default="uploader_config.toml",  # don't use
                        help="The path to the uploader config file "
                             "eg uploader_config.toml",
                        metavar="UPLOADER_CONFIG")

    return parser


def _set_default_options(options):
    """
    Set any options missing from the config file to their defaults.

    Options read from the config file are initially in the dictionary
    options.config. We move any 'commandline settable' key/value pairs
    from options.config to the options. There are also keys in options.config
    that can only be set via the config file - these aren't moved, but may be
    set to default values if they are missing and required.

    :param options: A dict-like object containing commandline and config file
                    options.
    :type options: dict | attrdict.AttrDict
    :return: The modified dict-like with defaults and overrides.
    :rtype: dict | attrdict.AttrDict
    """

    options_defaults = {
        'run_storage_base': None,
        'watch': False,
        'verbose': True,
        'uploader_config': 'uploader_config.toml',
        'logging_config': 'logging_config.toml',
        'notify_frequency': 60*24,  # daily
        'skip_bad_permissions': True,
    }

    for k, v in options_defaults.items():
        # Tranfer any known values set in options.config to the top level
        # options.
        # Any key not present in the config file gets set to the default value.
        if k not in options.config:
            options[k] = v
        else:
            options[k] = options.config[k]
            del options.config[k]

        if options[k] is None:
            options[k] = v

    return options


def _options_commandline_overrides(options):
    """
    Override any config file option with the value specified on the
    commandline.

    This is typically called after

    :param options: A dict-like object containing commandline and config file
                    options.
    :type options: dict | attrdict.AttrDict
    :return: The modified dict-like with defaults and overrides.
    :rtype: dict | attrdict.AttrDict
    """
    cmdline_values = {
        'run_storage_base': options.run_storage_base,
        'watch': options.watch,
        'verbose': options.verbose,
        'uploader_config': options.uploader_config,
        'logging_config': options.logging_config,
    }

    # Commandline options override any value in the config file.
    for k, v in cmdline_values.items():
        if v is not None:
            options[k] = v

    return options


def _attrdict_copy(d):
    """
    Copy an AttrDict (or create an AttrDict from a dict or argparse Namespace
    object.

    :param d:
    :type d:
    :return:
    :rtype:
    """
    if isinstance(d, collections.Iterable):
        return AttrDict(dict(d))
    else:
        return AttrDict(vars(d))


@atexit.register
def _set_current_task_status_on_exit():
    """
    Attempt to set current task status to 'error' if we have a premature
    exit.
    """
    global current
    global taskdb

    if current.task is not None \
            and taskdb is not None \
            and current.task.status != tasks.COMPLETE:
        current.task.status = tasks.ERROR
        taskdb.put(current.task.task_name, current.task)


def parse_run_specific_config(run_dir, options):
    """
    Override any options in options.config with values found in
    an autoprocess_config.toml file in the run directory.

    options.global_config is used as the set of initial default values.

    :param run_dir: The sequencing run directory.
    :type run_dir: str
    :param options: The options, containing a global_options attribute.
    :type options: AttrDict
    :return: A options object where options.config has value overridden by
             those found in the run directory config, if present.
    :rtype: AttrDict
    """
    run_config_path = path.join(run_dir, 'autoprocess_config.toml')
    if os.path.isfile(run_config_path):
        options.config = config_helper.get_config_toml(
            config_file=run_config_path,
            defaults=options.global_config)
        options.config = _attrdict_copy(options.config)
        logging.info("Applied run specific configuration from: %s",
                     run_config_path)
    else:
        options.config = _attrdict_copy(options.global_config)

    return options


def run_in_console():
    parser = setup_commandline_args()
    options = parser.parse_args()
    # options, argv = parser.parse_known_args()

    # TODO: Implement a support for global config (found in expected places, eg
    #       ~/.config/mytardis_ngs_ingestor/autoprocess_config.toml)
    #       which gets overridden by any options in a config found in the
    #       run directory. Might also allow multiple --config args to be used,
    #       with the second one listed taking precedence over the first
    #       (so we can compose a commandline with a global and run-local config
    #        and not require use of default locations)

    # Read the global autoprocessing config file, add values to
    # options.config dictionary
    options.config = AttrDict()
    if options.config_file and os.path.isfile(options.config_file):
        try:
            options.config = config_helper.get_config_toml(
                options.config_file,
                default_config_filename='autoprocess_config.toml')
            # options.config = _attrdict_copy(options.config)
            options = _attrdict_copy(options)
            options = _set_default_options(options)
            options = _options_commandline_overrides(options)
        except IOError as e:
            parser.error("Cannot read config file: %s" %
                         options.config_file)
            sys.exit(1)

    setup_logging(logging_config_file=options.logging_config)

    # If a config for illumina_uploader was provided, read it and assign it to
    # options.mytardis_uploader for use in the mytardis_upload task
    # options['config']['mytardis_uploader'] = None
    if options.uploader_config and os.path.isfile(options.uploader_config):
        uploader_only_parser = _add_uploader_config_argparser()
        uploader_config = _attrdict_copy(
            illumina_uploader.get_config_options(
                config_file_attr='uploader_config',
                ignore_unknown_args=True,
                parser=uploader_only_parser,
            )
        )
        options['config']['mytardis_uploader'] = uploader_config

    # options.global_config is a static copy of the config based
    # on the initial 'global' config file (default or --config)
    # and the commandline options. options.config may be
    # overridden by values in a run-specific config file, using
    # options.global_config as initial values.
    options.global_config = _attrdict_copy(options.config)

    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit()
    try:
        if options.run_storage_base and options.run_path:
            options.run_storage_base = None
            logging.warning("Please use only --runs or --single-run, not both.")

        if options.watch and not options.run_storage_base:
            logging.error("You must specify the --runs option if using --watch.")

        if options.run_storage_base:
            if options.watch:
                watch_runs(options.run_storage_base, options)
                sys.exit()
            else:
                logging.info('Running single autoprocess pass on: %s',
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
