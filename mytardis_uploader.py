# MyTardis Uploader v1.1
# Adopted and enhanced Andrew Perry <Andrew.Perry@monash.edu>
# Steve Androulakis <steve.androulakis@monash.edu>
# Thanks Grischa Meyer <grischa.meyer@monash.edu> for initial script

import logging
logger = logging.getLogger('mytardis_uploader')

import urllib2
import base64
import os
import mimetypes
import json
import requests
from requests.auth import HTTPBasicAuth
from time import strftime
import csv
from enum import Enum

class StorageMode(Enum):
    upload  =  1
    staging =  2
    shared  =  3

DEFAULT_STORAGE_MODE = StorageMode.upload

class PreemptiveBasicAuthHandler(urllib2.BaseHandler):
    def __init__(self, password_mgr=None):
        if password_mgr is None:
            password_mgr = urllib2.HTTPPasswordMgrWithDefaultRealm()
            self.passwd = password_mgr
            self.add_password = self.passwd.add_password

    def http_request(self, req):
        uri = req.get_full_url()
        user, pw = self.passwd.find_user_password(None, uri)
        if pw is None:
            return req

        raw = "%s:%s" % (user, pw)
        auth = 'Basic %s' % base64.b64encode(raw).strip()
        req.add_unredirected_header('Authorization', auth)
        return req


class MyTardisUploader:
    def __init__(self,
                 mytardis_url,
                 username,
                 password,
                 ):

        self.mytardis_url = mytardis_url
        self.v1_api_url = mytardis_url + "/api/v1/%s"
        self.username = username
        self.password = password

    def upload_directory(self,
                         file_path,
                         title='',
                         institute='',
                         description='',
                         test_run=False,
                         storage_mode=DEFAULT_STORAGE_MODE,
                         base_path=None,
                         exclude_patterns=[],
                         ):

        title = title or os.path.basename(os.path.abspath(file_path))

        logger.info('Creating experiment: %s', title)

        created = False
        exp_url = "/test/"
        if not test_run:
            agent_url = "https://github.com/steveandroulakis/mytardis-uploader"
            agent_footer = "Automatically generated on %s by %s" % \
                           (agent_url, strftime("%Y-%m-%d %H:%M:%S"))
            exp_url = self.create_experiment(title, institute or '', "%s\n\n%s"
                                             % (description or
                                                'No description.',
                                                agent_footer))
            created = True

        for item in os.listdir(file_path):
            full_path = os.path.join(file_path, item)

            if item.startswith('.') or item == 'metadata':
                continue  # filter files/dirs starting with .

            if any(regex.search(full_path) for regex in exclude_patterns):
                logger.info('Skipping excluded directory: %s', full_path)
                continue

            if os.path.isfile(full_path):
                logger.info('Skipping root-level file: %s', item)
            else:
                logger.info('Found directory: %s', item)

                logger.info('Creating dataset: %s', item)

                parameter_sets_list = \
                    self._get_dataset_parametersets_from_json(file_path,
                                                              item) or \
                    self._get_dataset_parametersets_from_csv(file_path, item)

                if parameter_sets_list:
                    logger.info("Found parameters for %s", item)

                dataset_url = "/__dry_run__/"
                if not test_run:
                    dataset_url = self.create_dataset(
                        '%s' % item,
                        [self._get_path_from_url(exp_url)],
                        parameter_sets_list
                    )

                for dirname, dirnames, filenames in os.walk(full_path):

                    for filename in filenames:
                        full_path = os.path.join(full_path, dirname, filename)

                        if filename.startswith('.'):
                            continue  # filter files/dirs starting with .

                        if any(regex.search(full_path)
                               for regex in exclude_patterns):
                            logger.info('Skipping excluded file: %s', full_path)
                            continue

                        logger.info("Uploading file: '%s' (dataset='%s')",
                                    full_path, item)

                        parameter_sets_list = \
                            self._get_datafile_parametersets_from_csv(full_path,
                                                                      item,
                                                                      filename)

                        if parameter_sets_list:
                            logger.info("Found parameters for %s", filename)

                        if not test_run:
                            # the replica_url should be the relative path
                            # to the file at the shared storage location
                            # (ie relative to the base path defined by
                            # StorageBox.location in the Django model)
                            # this applies to replica.protocol='file' locations
                            # but probably not S3/Swift style object store
                            # locations with real http urls.
                            replica_url = full_path
                            if storage_mode == StorageMode.shared:
                                if base_path is not None:
                                    # eg, if storage box base path is:
                                    # /data/bigstorage/
                                    # and absolute file path is
                                    # /data/bigstorage/expt1/dataset1/file.txt
                                    # then replica_url should be:
                                    # expt1/dataset1/file.txt
                                    replica_url = os.path.relpath(full_path,
                                                                  base_path)

                            self.upload_file(full_path,
                                             self._get_path_from_url(dataset_url),
                                             parameter_sets_list,
                                             storage_mode=storage_mode,
                                             storage_box_name='default',
                                             replica_url=replica_url)

        if created:
            exp_id = exp_url.rsplit('/')[-2]
            new_exp_url = "%s/experiment/view/%s/" % (self.mytardis_url, exp_id)
            logger.info("Experiment created: %s", new_exp_url)
            return new_exp_url

        else:
            logger.info("Dry run complete.")
            return "http://example.com/test/success"

    def _get_parametersets_from_csv(self,
                                    entity_type,
                                    file_path,
                                    filename):

        """
        entity_type can be 'datafile' or 'dataset'.
        """

        schema = "http://test.com"

        schema_path = '%s/metadata/schema.txt' % (os.path.abspath(file_path))

        try:
            with open(schema_path, 'rb') as schema_file:
                spamreader = csv.reader(schema_file,
                                        delimiter=',',
                                        quotechar='|')

                for row in spamreader:
                    if row[0].strip() == entity_type:
                        schema = row[1]
                        break
        except IOError:
            return []

        parameter_list = []

        try:
            with open(filename, 'rb') as csvfile:
                spamreader = csv.reader(csvfile,
                                        delimiter=',',
                                        quotechar='|')

                for row in spamreader:
                    parameter_list.append({u'name': row[0].strip(),
                                           u'value': row[1].strip()})
        except IOError:
            return []

        parameter_set = {'schema': schema,
                         'parameters': parameter_list}
        parameter_sets = [parameter_set]

        return parameter_sets

    def _get_dataset_parametersets_from_csv(self, file_path, dataset_name):
        filename = '%s/metadata/%s_metadata.csv' % (os.path.abspath(file_path),
                                                    dataset_name)

        parameter_sets = self._get_parametersets_from_csv('dataset',
                                                          file_path,
                                                          filename)

        return parameter_sets

    def _get_datafile_parametersets_from_csv(self,
                                             file_path,
                                             dataset_name,
                                             file_name):
        filename = '%s/metadata/%s_%s_metadata.csv' % \
                   (os.path.abspath(file_path),
                    dataset_name,
                    file_name)
        #print filename

        parameter_sets = self._get_parametersets_from_csv('datafile',
                                                          file_path,
                                                          filename)

        return parameter_sets

    # TODO _get_datafile_parametersets_from_json

    def _get_dataset_parametersets_from_json(self, file_path, dataset_name):
        filename = '%s/metadata/%s_metadata.json' % (os.path.abspath(file_path),
                                                     dataset_name)

        try:
            with open(filename) as f:
                parametersets = json.loads(f.read())
                return parametersets
        except IOError:
            pass

        return []

    def _do_create_request(self, data, urlend, method="POST"):
        url = self.v1_api_url % urlend
        headers = {'Accept': 'application/json',
                   'Content-Type': 'application/json'}

        auth_handler = PreemptiveBasicAuthHandler()
        auth_handler.add_password(realm=None,
                                  uri=url,
                                  user=self.username,
                                  passwd=self.password)
        opener = urllib2.build_opener(auth_handler)
        # ...and install it globally so it can be used with urlopen.
        urllib2.install_opener(opener)

        myrequest = urllib2.Request(url=url, data=data,
                                    headers=headers)
        myrequest.get_method = lambda: method
        # print myrequest.get_full_url() + " " + myrequest.data
        output = urllib2.urlopen(myrequest)

        return output

    def _md5_file_calc(self, file_path, blocksize=None):
        """
        Calculates the MD5 checksum of a file, returns the hex digest as a
        string. Streams the file in chunks of 'blocksize' to prevent running
        out of memory when working with large files.

        :param file_path: string
        :param blocksize: int
        :return: string
        """
        if not blocksize:
            blocksize = 128

        import hashlib
        md5 = hashlib.md5()
        with open(file_path, 'rb') as f:
            while True:
                chunk = f.read(blocksize)
                if not chunk:
                    break
                md5.update(chunk)
        return md5.hexdigest()

    def _send_datafile(self, data, urlend, filename=None):
        url = self.v1_api_url % urlend

        # we need to use requests_toolbelt here to prepare the multipart
        # encoded form data since vanilla requests can't stream files
        # when POSTing forms of this type and will run out of RAM
        # when encoding large files
        # See: https://github.com/kennethreitz/requests/issues/1584
        from requests_toolbelt import MultipartEncoder

        with open(filename, 'rb') as f:
            form = MultipartEncoder(fields={'json_data': data,
                                            'attached_file': (f, 'text/plain')})
            headers = {'Accept': 'application/json',
                       'Content-Type': form.content_type}

            response = requests.post(url,
                                     data=form,
                                     headers=headers,
                                     auth=HTTPBasicAuth(self.username,
                                                        self.password)
                                     )
            return response

    def _register_datafile_staging(self, data, urlend):
        raise NotImplementedError("Registering datafiles in a staging location"
                                  "is not currently implemented.")

    def _register_datafile_shared_storage(self, data, urlend):
        url = self.v1_api_url % urlend
        headers = {'Accept': 'application/json',
                   'Content-Type': 'application/json'}

        response = requests.post(url,
                                 data=data,
                                 headers=headers,
                                 auth=HTTPBasicAuth(self.username,
                                                    self.password)
                                 )
        return response


    def _get_header(self, headers, key):
        # from urllib2 style
        # get_header(headers, 'Location')

        import string

        location = None
        for header in string.split(headers, '\n'):
            if header.startswith('%s: ' % key):
                location = string.split(header, '%s: ' % key)[1].strip()
                break

        return location

    def _get_path_from_url(self, url_string):

        from urlparse import urlparse

        o = urlparse(url_string)

        return o.path

    def _format_parameter_set(self, schema, parameter_list):

        parameter_set = {'schema': schema, 'parameters': parameter_list}

        return parameter_set

    def create_experiment(self, title, institution, description,
                          author_list=None):

        # test authors ...
        # author_list = []
        # author_list.append({u'name': 'Daouda A.K. Traore', u'url': ''})
        # author_list.append({u'name': 'James C Whisstock', u'url': ''})

        if not author_list:
            author_list = []

        exp_dict = {
            u'description': description,
            u'institution_name': institution,
            u'title': title,
            u'authors': author_list
        }

        exp_json = json.dumps(exp_dict)

        # print exp_json

        data = self._do_create_request(exp_json, 'experiment/')

        return data.info().getheaders('Location')[0]

    def create_dataset(self, description, experiments_list,
                       parameter_sets_list=None, immutable=False):

        if not parameter_sets_list:
            parameter_sets_list = []

        dataset_dict = {
            u'description': description,
            u'experiments': experiments_list,
            u'immutable': immutable,
            u'parameter_sets': parameter_sets_list
        }

        dataset_json = json.dumps(dataset_dict)

        data = self._do_create_request(dataset_json, 'dataset/')

        return data.info().getheaders('Location')[0]

    def upload_file(self, file_path, dataset_path,
                    parameter_sets_list=None,
                    storage_mode=DEFAULT_STORAGE_MODE,
                    storage_box_name='default',
                    replica_url=''):
        # print upload_file('cli.py',
        #                   '/api/v1/dataset/143/').headers['location']

        if not parameter_sets_list:
            parameter_sets_list = []

        filename = os.path.basename(file_path)

        replica_list = [{u'url': replica_url,
                         u'location': storage_box_name,
                         u'protocol': u'file'},
                        ]
        file_size = os.path.getsize(file_path)
        # Hack to work around MyTardis not accepting
        # files of zero bytes
        #file_size = (file_size if file_size > 0 else -1)

        file_dict = {
            u'dataset': dataset_path,
            u'filename': filename,
            u'md5sum': self._md5_file_calc(file_path),
            u'mimetype': mimetypes.guess_type(file_path)[0],
            u'size': file_size,
            u'parameter_sets': parameter_sets_list,
            u'replicas': replica_list,
        }

        if storage_mode == StorageMode.shared:
            data = self._register_datafile_shared_storage(
                json.dumps(file_dict),
                'dataset_file/'
            )
        elif storage_mode == StorageMode.staging:
            data = self._register_datafile_staging(
                json.dumps(file_dict),
                'dataset_file/'
            )
        elif storage_mode == StorageMode.upload:
            delattr(file_dict, u'replicas')
            data = self._send_datafile(
                json.dumps(file_dict),
                'dataset_file/',
                filename=file_path
            )
        else:
            # we should never get here
            raise Exception("Invalid storage mode: " + storage_mode.name)

        location = getattr(data.headers, 'location', None)
        # print "Location: " + str(location)

        if (data.status_code > 499):
            logger.error("Registration of data file failed: %s", data.content)
            import sys
            sys.exit()

        return location


def run():
    ####
    # Le Script
    ####
    # steve.androulakis@monash.edu
    ####

    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.DEBUG)

    try:
        from colorlog import ColoredFormatter

        color_formatter = ColoredFormatter(
            '%(bg_white)s%(fg_black)s%(asctime)-8s%(reset)s\t'
            '%(log_color)s%(levelname)-8s%(reset)s\t'
            '%(white)s%(message)s',
            reset=True,
            log_colors={
                'DEBUG':    'cyan',
                'INFO':     'green',
                'WARNING':  'yellow',
                'ERROR':    'red',
                'CRITICAL': 'red,bg_white',
            },
            secondary_log_colors={},
            style='%'
        )
        console_handler.setFormatter(color_formatter)
    except:
        logging.basicConfig(format='%(asctime)s\t'
                                   # '%(name)-12s\t'
                                   '%(levelname)-8s\t'
                                   '%(message)s')
    logger.setLevel(logging.DEBUG)

    logger.addHandler(console_handler)
    logger.setLevel(logging.DEBUG)

    from argparse import ArgumentParser
    from appsettings import SettingsParser
    import getpass

    print """\
    MyTardis uploader generic v1
    Steve Androulakis <steve.androulakis@monash.edu>
    Uploads the given directory as an Experiment, and the immediate
    sub-directories below it as Datasets in MyTardis.

    eg. python mytardis_uploader.py -l http://mytardis-server.com.au -u steve \
                                    -f /Users/steve/Experiment1/"

    If present, metadata is harvested from:
      <path>/metadata/schema.txt
      <path>/metadata/<dataset_name>_metadata.csv
      <path>/metadata/<dataset_name>_metadata.json
      <path>/metadata/<dataset_name>_<filename>_metadata.csv
      <path>/metadata/<dataset_name>_<filename>_metadata.json

    """

    def add_args(parser):
        """
        Takes an argparse.ArgumentParser-like object and adds commandline
        parameters to detect and capture values from.
        :param parser: ArgumentParser
        :return:
        """
        parser.add_argument('--config',
                            dest="config_file",
                            type=str,
                            metavar="MYTARDIS_UPLOADER_CONFIG")
        parser.add_argument("-f", "--path",
                            dest="path",
                            type=str,
                            help="The PATH of the experiment to be uploaded",
                            metavar="PATH")
        parser.add_argument("--storage-base-path",
                            dest="storage_base_path",
                            type=str,
                            help="The STORAGE_BASE_PATH of all experiments,"
                                 "when using 'shared' storage mode.",
                            metavar="STORAGE_BASE_PATH")
        parser.add_argument("-l", "--url",
                            dest="url",
                            type=str,
                            help="The URL to the MyTardis installation",
                            metavar="URL")
        parser.add_argument("-u", "--username",
                            dest="username",
                            type=str,
                            help="Your MyTardis USERNAME",
                            metavar="USERNAME")
        parser.add_argument("--password",
                            dest="password",
                            type=str,
                            help="You should probably never use this option from "
                                 "the command line. Be sensible.",
                            metavar="PASSWORD")
        parser.add_argument("-t", "--title",
                            dest="title",
                            type=str,
                            help="Experiment TITLE",
                            metavar="TITLE")
        parser.add_argument("-d", "--description",
                            dest="description",
                            type=str,
                            help="Experiment DESCRIPTION",
                            metavar="DESCRIPTION")
        parser.add_argument("-i", "--institute",
                            dest="institute",
                            type=str,
                            help="Experiment INSTITUTE (eg university)",
                            metavar="INSTITUTE")
        parser.add_argument("-r", "--dry",
                            action="store_true",
                            dest="dry_run",
                            default=False,
                            help="Dry run (don't create anything)")
        parser.add_argument("--storage-mode",
                            dest="storage_mode",
                            type=str,
                            metavar="STORAGE_MODE",
                            default='upload',
                            help="Specify if the data files are to be uploaded, "
                                 "or registered in the database at a staging or "
                                 "shared storage area without uploading. "
                                 "Valid values are: upload, staging or shared."
                                 "Defaults to upload.")
        parser.add_argument("--exclude",
                            dest="exclude",
                            action="append",
                            help="Exclude files with paths matching this regex. "
                                 "Can be specified multiple times.",
                            metavar="REGEX")


    def get_config(default_config_filename='uploader_config.yaml'):
        """
        Parses a config file (default or commandline specified), then
        overrides any settings with those specified on the command line.
        Returns the appsetting.SettingsParser instance and the config options
        object (argparse-like).

        We do something a little unusual here to allow a --config option
        while using the appsettings module.
        First we parse the commandline using standard argparse and look for
        the --config option. If we find it, we read the config via
        appsettings.SettingsParser. Then we reparse the commandline using that
        SettingsParser to override any config file options with commandline
        options specified.

        :return: parser, options
        """
        preparser = ArgumentParser()
        add_args(preparser)
        precheck_options = preparser.parse_args()

        parser = None

        if precheck_options.config_file:
            try:
                with open(precheck_options.config_file, 'r') as f:
                    parser = SettingsParser(yaml_file=f)
            except IOError:
                preparser.error("Cannot read config file: %s" %
                                precheck_options.config_file)
        elif os.path.isfile(default_config_filename):
            with open(default_config_filename, 'r') as f:
                parser = SettingsParser(yaml_file=f)
        else:
            parser = SettingsParser()

        add_args(parser)
        options = parser.parse_args()

        return parser, options

    parser, options = get_config()

    if not options.path:
        parser.error('File path not given')

    if not options.url:
        parser.error('URL to MyTardis instance not given')

    if not options.username:
        parser.error('MyTardis username not given')

    # valid_storage_modes = ['upload', 'staging', 'shared']
    valid_storage_modes = [m.name for m in StorageMode]
    if options.storage_mode and \
       options.storage_mode not in valid_storage_modes:
        parser.error('--storage-mode must be one of: ' +
                     ', '.join(valid_storage_modes))

    if options.storage_mode is 'shared' and not options.storage_base_path:
        parser.error("--storage-base-path (storage_base_path) must be"
                     "specified when using 'shared' storage mode.")

    exclude_patterns = []
    if options.exclude:
        import re

        for regex in options.exclude:
            # strip matching quotes around regex if present
            if regex.startswith('"') and regex.endswith('"'):
                regex = regex[1:-1]
            elif regex.startswith("'") and regex.endswith("'"):
                regex = regex[1:-1]
            exclude_patterns.append(re.compile(regex))

        logger.info("Ignoring files that match: %s\n",
                    ' | '.join(options.exclude))

    pw = options.password
    if not pw:
        pw = getpass.getpass()

    file_path = options.path
    if file_path is '.':
        file_path = os.getcwd()
    title = options.title
    institute = options.institute
    description = options.description
    test_run = options.dry_run
    mytardis_url = options.url
    username = options.username
    password = pw
    storage_mode = StorageMode[options.storage_mode]

    mytardis_uploader = MyTardisUploader(mytardis_url,
                                         username,
                                         password)

    mytardis_uploader.upload_directory(file_path,
                                       title=title,
                                       description=description,
                                       institute=institute,
                                       test_run=test_run,
                                       storage_mode=storage_mode,
                                       base_path=options.storage_base_path,
                                       exclude_patterns=exclude_patterns)


if __name__ == "__main__":
    run()
