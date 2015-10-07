#!/usr/bin/env python
# MyTardis Uploader
# Adopted and enhanced Andrew Perry <Andrew.Perry@monash.edu>
# Steve Androulakis <steve.androulakis@monash.edu>
# Thanks Grischa Meyer <grischa.meyer@monash.edu> for initial script

from __init__ import __version__
import logging

logger = logging.getLogger('mytardis_uploader')

import os
import sys
import mimetypes
import json
import requests
from requests.auth import HTTPBasicAuth
import backoff
from time import strftime
import datetime
import csv

import urllib3
logging.captureWarnings(True)

# https://urllib3.readthedocs.org/en/latest/contrib.html#module-urllib3.contrib.pyopenssl
try:
    import urllib3.contrib.pyopenssl
    urllib3.contrib.pyopenssl.inject_into_urllib3()
except ImportError:
    pass

DEFAULT_STORAGE_MODE = 'upload'


class MyTardisUploader:
    def __init__(self,
                 mytardis_url,
                 username,
                 password,
                 storage_mode=DEFAULT_STORAGE_MODE,
                 storage_box_location='',
                 storage_box_name='default',
                 verify_certificate=True,
                 ):

        self.mytardis_url = mytardis_url
        self.v1_api_url = mytardis_url + "/api/v1/%s"
        self.username = username
        self.password = password
        self.storage_mode = storage_mode
        self.storage_box_location = storage_box_location
        self.storage_box_name = storage_box_name
        self.user_agent_url = "https://github.com/pansapiens/mytardis-uploader"
        self.user_agent_name = "MyTardisUploader"
        self.user_agent = "%s/%s (%s)" % (self.user_agent_name,
                                          __version__,
                                          self.user_agent_url)
        # True, False, or the path to the certificate (.pem)
        self.verify_certificate = verify_certificate

    def _json_request_headers(self):
        return {'Accept': 'application/json',
                'Content-Type': 'application/json',
                'User-Agent': self.user_agent}

    @staticmethod
    def dict_to_json(d):
        """
        Serialize a dictionary to JSON, correctly handling datetime.datetime
        objects (to ISO 8601 dates, as strings).

        :type d: dict
        :rtype: str
        """
        date_handler = lambda obj: (
            obj.isoformat(' ')
            if isinstance(obj, datetime.date)
            or isinstance(obj, datetime.datetime)
            else None
        )
        return json.dumps(d, default=date_handler)

    def upload_directory(self,
                         file_path,
                         title='',
                         institute='',
                         description='',
                         test_run=False,
                         exclude_patterns=None,
                         ):

        title = title or os.path.basename(os.path.abspath(file_path))

        logger.info('Creating experiment: %s', title)

        created = False
        exp_url = "/test/"
        if not test_run:
            # "2015-05-19 00:52:38.612843"
            # start_time = datetime.datetime.now().isoformat(' ')
            # end_time = datetime.datetime.now().isoformat(' ')

            agent_footer = "Automatically generated by %s on %s" % \
                           (self.user_agent, strftime("%Y-%m-%d %H:%M:%S"))
            expt_dict = {'title': title,
                         'institute_name': institute or '',
                         'description': "%s\n\n%s" % (description or
                                                      'No description.',
                                                      agent_footer),
                         'start_time': datetime.datetime.now(),
                         'end_time': datetime.datetime.now(),
            }
            exp_url = self.create_experiment(expt_dict)
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
                    dataset_dict = {
                        'description': '%s' % item,
                        'experiments': [self._get_path_from_url(exp_url)],
                        'parameter_sets': parameter_sets_list
                    }
                    dataset_url = self.create_dataset(dataset_dict)

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
                            if self.storage_mode == 'shared':
                                if self.storage_box_location:
                                    # eg, if storage box base path is:
                                    # /data/bigstorage/
                                    # and absolute file path is
                                    # /data/bigstorage/expt1/dataset1/file.txt
                                    # then replica_url should be:
                                    # expt1/dataset1/file.txt
                                    replica_url = os.path.relpath(
                                        full_path,
                                        self.storage_box_location
                                    )

                            self.upload_file(full_path,
                                             self._get_path_from_url(
                                                 dataset_url
                                             ),
                                             parameter_sets_list,
                                             replica_url=replica_url)

        if created:
            exp_id = exp_url.rsplit('/')[-2]
            new_exp_url = "%s/experiment/view/%s/" % (self.mytardis_url, exp_id)
            logger.info("Experiment created: %s", new_exp_url)
            return new_exp_url

        else:
            logger.info("Dry run complete.")
            return "http://example.com/test/success"

    def _get_parametersets_from_csv(self, entity_type,
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
        # print filename

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

    def _raise_502(self, response):
        e = requests.exceptions.RequestException(response=response)
        e.message = "%s %s" % (response.status_code, response.reason)
        raise e

    @backoff.on_exception(backoff.expo,
                          requests.exceptions.RequestException,
                          max_tries=8)
    def _do_post_request(self, data, urlend):
        url = self.v1_api_url % (urlend.rstrip('/') + '/')
        headers = self._json_request_headers()

        try:
            response = requests.post(url,
                                     data=data,
                                     headers=headers,
                                     auth=HTTPBasicAuth(self.username,
                                                        self.password),
                                     verify=self.verify_certificate,
                                     )
            # 502 Bad Gateway triggers retries, since the proxy web
            # server (eg Nginx or Apache) in front of MyTardis could be
            # temporarily restarting
            if response.status_code == 502:
                self._raise_502(response)

        except requests.exceptions.RequestException, e:
            logger.error("Request failed : %s : %s", e.message, url)
            raise e

        return response

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

    @backoff.on_exception(backoff.expo,
                          requests.exceptions.RequestException,
                          max_tries=8)
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
                                            'attached_file': ('text/plain', f)})
            headers = self._json_request_headers()
            headers['Content-Type'] = form.content_type

            try:
                response = requests.post(url,
                                         data=form,
                                         headers=headers,
                                         auth=HTTPBasicAuth(self.username,
                                                            self.password),
                                         verify=self.verify_certificate,
                                         )

                if response.status_code == 502:
                    self._raise_502(response)

            except requests.exceptions.RequestException, e:
                logger.error("Request failed : %s : %s", e.message, url)
                raise e

            return response

    @backoff.on_exception(backoff.expo,
                          requests.exceptions.RequestException,
                          max_tries=8)
    def _register_datafile_staging(self, data, urlend):
        raise NotImplementedError("Registering datafiles in a staging location"
                                  "is not currently implemented.")

    @backoff.on_exception(backoff.expo,
                          requests.exceptions.RequestException,
                          max_tries=8)
    def _register_datafile_shared_storage(self, data, urlend):
        url = self.v1_api_url % urlend
        headers = self._json_request_headers()

        try:
            response = requests.post(url,
                                     data=data,
                                     headers=headers,
                                     auth=HTTPBasicAuth(self.username,
                                                        self.password),
                                     verify=self.verify_certificate,
                                     )

            if response.status_code == 502:
                self._raise_502(response)

        except requests.exceptions.RequestException, e:
            logger.error("Request failed : %s : %s", e.message, url)
            raise e

        return response

    def _get_path_from_url(self, url_string):
        from urlparse import urlparse
        o = urlparse(url_string)
        return o.path

    def _format_parameter_set(self, schema, parameter_list):
        parameter_set = {'schema': schema, 'parameters': parameter_list}
        return parameter_set

    def create_experiment(self, expt_dict):
        """
        Execute an HTTP request to create an experiment on the server.

        :param expt_dict: A dictionary representing the experiment, ready to be
                          serialized to JSON for the MyTardis REST API.
        :type expt_dict: dict
        :return: The url path of the created experiment.
        :rtype: str
        """

        if not expt_dict.get('created_time', None):
            from datetime import datetime
            # Z is UTC
            expt_dict['created_time'] = datetime.utcnow()  # .isoformat('Z')

        # author_list is of the format:
        # [{u'name': 'Daouda A.K. Traore', u'url': ''},
        #  {u'name': 'James C Whisstock', u'url': ''}]
        # if author_list and isinstance(author_list, list):
        #     exp_dict[u'authors'] = author_list

        expt_json = self.dict_to_json(expt_dict)

        # logger.debug("Experiment JSON: %s", expt_json)

        data = self._do_post_request(expt_json, 'experiment')

        if not data.ok or 'Location' not in data.headers:
            logger.error("Creating experiment failed: %s (%s) - %s",
                         data.reason,
                         data.status_code,
                         data.content)
            sys.exit(1)

        return data.headers.get('Location', None)

    def query_instrument(self, name):
        url = self.v1_api_url % 'instrument'
        query_params = {u'name': name}
        headers = self._json_request_headers()
        response = requests.get(url,
                                params=query_params,
                                headers=headers,
                                auth=HTTPBasicAuth(self.username,
                                                   self.password),
                                verify=self.verify_certificate,
                                )
        return response.json()

    def query_group(self, name):
        url = self.v1_api_url % 'group'
        query_params = {u'name': name}
        headers = self._json_request_headers()
        response = requests.get(url,
                                params=query_params,
                                headers=headers,
                                auth=HTTPBasicAuth(self.username,
                                                   self.password),
                                verify=self.verify_certificate,
                                )
        return response.json()

    def create_dataset(self, dataset, instrument=None):
        """
        Execute an HTTP request to create a Dataset on the server.

        :param dataset: A dictionary representing the dataset ready to be
                          serialized to JSON for the MyTardis REST API.
        :type dataset: dict
        :param instrument: An instrument name known by the server
        :type instrument: unicode
        :return: The url path of the dataset created
        :rtype: str
        """
        if instrument is not None:
            instrument_resource_uri = \
                self.query_instrument(instrument)['objects'][0]['resource_uri']
            dataset['instrument'] = instrument_resource_uri

        dataset_json = MyTardisUploader.dict_to_json(dataset)
        data = self._do_post_request(dataset_json, 'dataset')

        if not data.ok or 'Location' not in data.headers:
            logger.error("Creating dataset failed: %s", data.text)
            sys.exit(1)

        return data.headers.get('Location', None)

    def upload_file(self, file_path, dataset_url_path,
                    parameter_sets_list=None,
                    replica_url='',
                    md5_checksum=None):

        if not parameter_sets_list:
            parameter_sets_list = []

        filename = os.path.basename(file_path)

        replica_list = [{u'url': replica_url,
                         u'location': self.storage_box_name,
                         u'protocol': u'file'},
                        ]
        file_size = os.path.getsize(file_path)
        # Hack to work around MyTardis not accepting
        # files of zero bytes
        # file_size = (file_size if file_size > 0 else -1)
        if not md5_checksum:
            md5_checksum = self._md5_file_calc(file_path)

        file_dict = {
            u'dataset': dataset_url_path,
            u'filename': filename,
            u'md5sum': md5_checksum,
            u'mimetype': mimetypes.guess_type(file_path)[0],
            u'size': file_size,
            u'parameter_sets': parameter_sets_list,
            u'replicas': replica_list,
        }

        if self.storage_mode == 'shared':
            data = self._register_datafile_shared_storage(
                self.dict_to_json(file_dict),
                'dataset_file/'
            )
        elif self.storage_mode == 'staging':
            data = self._register_datafile_staging(
                self.dict_to_json(file_dict),
                'dataset_file/'
            )
        elif self.storage_mode == 'upload':
            # file_dict.pop(u'replicas', None)
            data = self._send_datafile(
                self.dict_to_json(file_dict),
                'dataset_file/',
                filename=file_path
            )
        else:
            # we should never get here
            raise Exception("Invalid storage mode: " + self.storage_mode)

        if not data.ok or 'Location' not in data.headers:
            logger.error("Registration of data file failed: %s", data.text)
            sys.exit(1)

        return data.headers.get('Location', None)

    def _resource_uri_to_id(self, uri):
        """
        Takes resource URI like: http://example.org/api/v1/experiment/998
        and returns just the id value (998).

        :type uri: str
        :rtype: int
        """
        from urlparse import urlparse
        resource_id = int(urlparse(uri).path.rstrip('/').split('/')[-1:][0])
        return resource_id

    def share_experiment_with_group(self, experiment, group_name):
        """
        Executes an HTTP request to share an experiment with a group,
        via updating the ObjectACL.

        :param experiment:  The url path to the experiment to update.
        :param group_name: Name of the group we with share with.
        :type experiment: str
        :type group_name: str
        :return: A requests Response object
        :rtype: Response
        """
        group_id = self.query_group(group_name)['objects'][0]['id']
        experiment_id = self._resource_uri_to_id(experiment)
        data = {
            u'pluginId': u'django_group',
            u'entityId': unicode(group_id),
            u'content_object': unicode(experiment),
            u'content_type': u'experiment',
            u'object_id': unicode(experiment_id),
            u'aclOwnershipType': 1,
            u'isOwner': True,
            u'canRead': True,
            u'canWrite': True,
            u'canDelete': False,
            u'effectiveDate': None,
            u'expiryDate': None
        }

        response = self._do_post_request(self.dict_to_json(data), 'objectacl')
        return response


def setup_logging(loglevel=logging.DEBUG):
    logger.setLevel(loglevel)

    console_handler = logging.StreamHandler(sys.stderr)
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
    except ImportError:
        logging.basicConfig(format='%(asctime)s\t'
                                   # '%(name)-12s\t'
                                   '%(levelname)-8s\t'
                                   '%(message)s')

    logger.addHandler(console_handler)

    low_level_request_logging = False
    if low_level_request_logging:
        try:
            import http.client as http_client
        except ImportError:
            # Python 2
            import httplib as http_client

        http_client.HTTPConnection.debuglevel = 1

        requests_log = logging.getLogger("requests.packages.urllib3")
        requests_log.setLevel(logging.DEBUG)
        requests_log.propagate = True

    return logger


def add_config_args(parser, add_extra_options_fn=None):
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
    parser.add_argument("--storage-box-name",
                        dest="storage_box_name",
                        type=str,
                        default='default',
                        help="The name of the MyTardis StorageBox to ingest"
                             "data into.",
                        metavar="STORAGE_BOX_NAME")
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
    parser.add_argument("--owners",
                        dest="experiment_owner_groups",
                        nargs="+",
                        help="A list of groups (by name) that will be given "
                             "ownership of ingested experiments (space"
                             "separated)",
                        metavar="EXPERIMENT_OWNER_GROUPS")
    parser.add_argument("-r", "--dry",
                        action="store_true",
                        dest="dry_run",
                        default=False,
                        help="Dry run (don't create anything)")
    parser.add_argument("--fast",
                        action="store_true",
                        dest="fast",
                        default=False,
                        help="Skip some time consuming steps but upload "
                             "incomplete metadata (eg, no md5 checksums)")
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
    # NOTE: When using this option in code, you probably want to use the value
    #       of options.verify_certificate rather than options.certificate
    #       (see logic in validate_config where the verify_certificate
    #        attribute is added)
    parser.add_argument("--certificate",
                        dest="certificate",
                        type=str,
                        default='True',
                        help="The SSL/TLS certificate used by the MyTardis "
                             "server. Required to estabilish the identity of "
                             "the server when using HTTPS and a self-signed "
                             "certificate."
                             "Valid values are a path to a certificate (.pem) "
                             "or INSECURE (all-caps). "
                             "If set to INSECURE, self-signed certificates are "
                             "blindly accepted. The identity of the server "
                             "cannot be guaranteed, increasing the likelyhood "
                             "of MITM attacks.",
                        metavar="CERTIFICATE")

    if add_extra_options_fn:
        add_extra_options_fn(parser)


def get_config(default_config_filename='uploader_config.yaml',
               add_extra_options_fn=None):
    """
    Parses a config file (default or commandline specified), then
    overrides any settings with those specified on the command line.
    Returns the appsetting.SettingsParser instance and the config options
    object (argparse.Namespace).

    We do something a little unusual here to allow a --config option
    while using the appsettings module.
    First we parse the commandline using standard argparse and look for
    the --config option. If we find it, we read the config via
    appsettings.SettingsParser. Then we reparse the commandline using that
    SettingsParser to override any config file options with commandline
    options specified.

    :param default_config_filename: str
    :return:
    :rtype: (argparse.ArgumentParser, object)
    """
    from argparse import ArgumentParser
    from appsettings import SettingsParser

    preparser = ArgumentParser()
    add_config_args(preparser,
                    add_extra_options_fn=add_extra_options_fn)
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

    add_config_args(parser,
                    add_extra_options_fn=add_extra_options_fn)

    options = parser.parse_args()

    return parser, options


def validate_config(parser, options):
    """
    Validates config options, throws errors and stops excution if
    there is an issue (invalid value or required value missing).

    :rtype : object
    :param options: object
    :param parser: argparse.ArgumentParser
    :return:
    """
    if not options.path:
        parser.error('File path not given')

    if not options.url:
        parser.error('URL to MyTardis instance not given')

    if not options.username:
        parser.error('MyTardis username not given')

    valid_storage_modes = ['upload', 'staging', 'shared']
    if options.storage_mode and \
       options.storage_mode not in valid_storage_modes:
        parser.error('--storage-mode must be one of: ' +
                     ', '.join(valid_storage_modes))

    if options.storage_mode == 'shared' and not options.storage_base_path:
        parser.error("--storage-base-path (storage_base_path) must be"
                     "specified when using 'shared' storage mode.")

    # We want to force certificate verification if this value is unset.
    # We set options.verify_certificate (a bool OR str) based on the
    # value of options.certificate.
    # An empty value in the YAML config gets serialized to the string 'None'
    options.verify_certificate = options.certificate
    if options.certificate == 'None' or \
       options.certificate.strip() == '':
        options.verify_certificate = True
    if options.certificate == 'INSECURE':
        options.verify_certificate = False


def get_exclude_patterns_as_regex_list(exclude_patterns=None):
    """
    Takes a list of strings are returns a list of compiled regexes.
    Strips enclosing quotes if present.

    :param options: list[str]
    :return:
    :rtype: list[re.__Regex]
    """
    exclude_regexes = []
    if exclude_patterns:
        import re

        for regex in exclude_patterns:
            # strip matching quotes around regex if present
            if regex.startswith('"') and regex.endswith('"'):
                regex = regex[1:-1]
            elif regex.startswith("'") and regex.endswith("'"):
                regex = regex[1:-1]
            exclude_regexes.append(re.compile(regex))

        logger.info("Ignoring files that match: %s\n",
                    ' | '.join(exclude_patterns))

    return exclude_regexes


def run():
    ####
    # Le Script
    ####
    # steve.androulakis@monash.edu
    ####
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

    setup_logging()

    parser, options = get_config()

    validate_config(parser, options)

    exclude_patterns = \
        get_exclude_patterns_as_regex_list(options.exclude)

    password = options.password
    if not password:
        password = getpass.getpass()

    file_path = options.path
    if file_path is '.':
        file_path = os.getcwd()
    title = options.title
    institute = options.institute
    description = options.description
    test_run = options.dry_run

    mytardis_uploader = MyTardisUploader(
        options.url,
        options.username,
        password,
        storage_mode=options.storage_mode,
        storage_box_location=options.storage_base_path,
        storage_box_name=options.storage_box_name,
        verify_certificate=options.verify_certificate)

    mytardis_uploader.upload_directory(
        file_path,
        title=title,
        description=description,
        institute=institute,
        test_run=test_run,
        exclude_patterns=exclude_patterns)

    # raw_data_expt = create_experiment
    # upload_directory_as_child_experiments(raw_data_expt)
    #   create_project_fastq_expt

if __name__ == "__main__":
    run()
