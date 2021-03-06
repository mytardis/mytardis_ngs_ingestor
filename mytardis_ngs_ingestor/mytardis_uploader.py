#!/usr/bin/env python
# MyTardis Uploader
# Adopted and enhanced Andrew Perry <Andrew.Perry@monash.edu>
# Steve Androulakis <steve.androulakis@monash.edu>
# Thanks Grischa Meyer <grischa.meyer@monash.edu> for initial script

from __future__ import print_function, absolute_import, division

__author__ = 'Andrew Perry <Andrew.Perry@monash.edu.au>'

from builtins import (bytes, str, open, super, range,
                      zip, round, input, int, pow, object)

try:
    from __init__ import __version__
except ImportError:
    from mytardis_ngs_ingestor.__init__ import __version__

import logging

import six
from six.moves.urllib.parse import urlparse, urljoin
import os
import sys
import mimetypes
import json
import requests
from requests.auth import AuthBase, HTTPBasicAuth
import backoff
from time import strftime
import datetime
import csv
import subprocess
import hashlib
from argparse import ArgumentParser
from appsettings import SettingsParser
from toolz.dicttoolz import merge as merge_dicts

from mytardis_ngs_ingestor.utils import config_helper, is_in_tree

import urllib3
logging.captureWarnings(True)

# https://urllib3.readthedocs.org/en/latest/contrib.html#module-urllib3.contrib.pyopenssl
try:
    import urllib3.contrib.pyopenssl
    urllib3.contrib.pyopenssl.inject_into_urllib3()
except ImportError:
    pass

DEFAULT_STORAGE_MODE = 'upload'


class TastyPieAuth(AuthBase):
    """
    Attaches HTTP headers for Tastypie API key Authentication to the given
    Request object.
    """
    def __init__(self, username, api_key):
        self.username = username
        self.api_key = api_key

    def __call__(self, r):
        r.headers['Authorization'] = 'ApiKey %s:%s' % (self.username,
                                                       self.api_key)
        return r


class MyTardisUploader:
    user_agent_name = __name__
    user_agent_url = 'https://github.com/mytardis/mytardis_ngs_ingestor'

    def __init__(self,
                 mytardis_url,
                 username,
                 password=None,
                 api_key=None,
                 storage_mode=DEFAULT_STORAGE_MODE,
                 storage_box_location='',
                 storage_box_name='default',
                 verify_certificate=True,
                 fast_mode=False,
                 exclude_patterns=None,
                 ):

        self.mytardis_url = mytardis_url
        # trailing slash matters
        self.v1_api_url = urljoin(mytardis_url, '/api/v1/%s/')
        self.username = username
        self.password = password
        self.api_key = api_key
        self.storage_mode = storage_mode
        self.storage_box_location = storage_box_location
        self.storage_box_name = storage_box_name
        self.user_agent = '%s/%s (%s)' % (MyTardisUploader.user_agent_name,
                                          __version__,
                                          MyTardisUploader.user_agent_url)
        # True, False, or the path to the certificate (.pem)
        self.verify_certificate = verify_certificate
        self.fast_mode = fast_mode
        self.exclude_patterns = \
            get_exclude_patterns_as_regex_list(exclude_patterns)
        if exclude_patterns is None:
            self.exclude_pattens = []

        if self.api_key is not None:
            self.auth = TastyPieAuth(self.username, self.api_key)
        elif self.password is not None:
            self.auth = HTTPBasicAuth(self.username, self.password)
        else:
            self.auth = None

    def __call__(self, file_path, dataset_url_path, **kwargs):
        """
        Allows an instance of MyTardisUploader to be called as a function
        that uploads a file. This primarily exists to allow
        the MyTardisUploader.upload_file method to be passed to
        `multiprocessing` (and `concurrent.futures.ProcessPoolExecutor`),
        since bound methods aren't pickleable.

        eg:

        uploader = MyTardisUploader()

        datafile_url = uploader("some_file.txt", "/dataset/99")

        An alternative solution that could be used with `multiprocessing` is
        to pickle the upload_file method like this: https://goo.gl/4uYVqN

        :param file_path:
        :type file_path: str
        :param dataset_url_path:
        :type dataset_url_path: str
        :param kwargs:
        :type kwargs:
        :return: A tuple of the filepath, dataset URL and datafile URL.
        :rtype: tuple(str, str, str)
        """
        return (file_path,
                dataset_url_path,
                self.upload_file(file_path, dataset_url_path, **kwargs))

    def _json_request_headers(self):
        return {'Accept': 'application/json',
                'Content-Type': 'application/json',
                'User-Agent': self.user_agent}

    # TODO: Replace use of this with the jsondate package.
    #       jsondate is now in requirements.txt
    #       use: import jsondate as json
    #       then just use json.dumps as normal, remove this method
    @staticmethod
    def dict_to_json(d):
        """
        Serialize a dictionary to JSON, correctly handling datetime.datetime
        objects (to ISO 8601 dates, as strings).

        :type d: dict
        :rtype: str
        """
        if not isinstance(d, dict):
            raise TypeError("Must be a dictionary")

        date_handler = lambda obj: (
            obj.isoformat(' ')
            if isinstance(obj, datetime.date)
            or isinstance(obj, datetime.datetime)
            else None
        )
        return json.dumps(d, default=date_handler)

    def upload_directory(self,
                         base_path,
                         title='',
                         institute='',
                         description='',
                         test_run=False):

        title = title or os.path.basename(os.path.abspath(base_path))

        logging.info('Creating experiment: %s', title)

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

        for item in os.listdir(base_path):
            full_path = os.path.join(base_path, item)

            if item.startswith('.') or item == 'metadata':
                continue  # filter files/dirs starting with .

            if any(regex.search(full_path) for regex in self.exclude_patterns):
                logging.info('Skipping excluded directory: %s', full_path)
                continue

            if os.path.isfile(full_path):
                logging.info('Skipping root-level file: %s', item)
            else:
                logging.info('Found directory: %s', item)

                logging.info('Creating dataset: %s', item)

                parameter_sets_list = \
                    self._get_dataset_parametersets_from_json(base_path,
                                                              item) or \
                    self._get_dataset_parametersets_from_csv(base_path, item)

                if parameter_sets_list:
                    logging.info("Found parameters for %s", item)

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
                               for regex in self.exclude_patterns):
                            logging.info('Skipping excluded file: %s', full_path)
                            continue

                        logging.info("Uploading file: '%s' (dataset='%s')",
                                    full_path, item)

                        parameter_sets_list = \
                            self._get_datafile_parametersets_from_csv(full_path,
                                                                      item,
                                                                      filename)

                        if parameter_sets_list:
                            logging.info("Found parameters for %s", filename)

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
            logging.info("Experiment created: %s", new_exp_url)
            return new_exp_url

        else:
            logging.info("Dry run complete.")
            return "http://example.com/test/success"

    def _get_parametersets_from_csv(self, entity_type, file_path, filename):

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

    def _get_parametersets_from_json(self, filename):
        try:
            with open(filename) as f:
                parametersets = json.loads(f.read())
                return parametersets
        except IOError:
            pass

        return []

    def _get_dataset_parametersets_from_json(self, file_path, dataset_name):
        filename = '%s/metadata/%s_metadata.json' % (os.path.abspath(file_path),
                                                     dataset_name)

        return self._get_parametersets_from_json(filename)

    def _get_datafile_parametersets_from_json(self,
                                              file_path,
                                              dataset_name,
                                              file_name):
        filename = '%s/metadata/%s_%s_metadata.json' % (
            os.path.abspath(file_path),
            dataset_name,
            file_name
        )

        return self._get_parametersets_from_json(filename)

    def _raise_request_exception(self, response):
        e = requests.exceptions.RequestException(response=response)
        e.message = "%s %s" % (response.status_code, response.reason)
        raise e

    # @backoff.on_exception(backoff.expo,
    #                       requests.exceptions.RequestException,
    #                       max_tries=8)
    # def do_get_request(self, action, query_params, extra_headers={}):
    #     url = self.v1_api_url % action
    #     headers = self._json_request_headers()
    #     headers = merge_dicts(headers, extra_headers)
    #
    #     try:
    #         response = requests.get(url,
    #                                 params=query_params,
    #                                 headers=headers,
    #                                 auth=self.auth,
    #                                 verify=self.verify_certificate,
    #                                 )
    #
    #         # 502 Bad Gateway triggers retries, since the proxy web
    #         # server (eg Nginx or Apache) in front of MyTardis could be
    #         # temporarily restarting
    #         if response.status_code == 502:
    #             self._raise_502(response)
    #
    #     except requests.exceptions.RequestException, e:
    #         logging.error("Request failed : %s : %s", e.message, url)
    #         raise e
    #
    #     return response

    def do_get_request(self, action, params, extra_headers=None):
        return self._do_request('GET', action,
                                params=params,
                                extra_headers=extra_headers)

    def do_post_request(self, action, data, extra_headers=None):
        return self._do_request('POST', action,
                                data=data,
                                extra_headers=extra_headers)

    @backoff.on_exception(backoff.expo,
                          requests.exceptions.RequestException,
                          max_tries=8)
    def _do_request(self, method, action,
                    data=None, params=None,
                    extra_headers=None,
                    api_url_template=None):

        if api_url_template is None:
            api_url_template = self.v1_api_url

        url = api_url_template % action

        headers = self._json_request_headers()
        if extra_headers is not None:
            headers = merge_dicts(headers, extra_headers)

        try:
            response = requests.request(method,
                                        url,
                                        data=data,
                                        params=params,
                                        headers=headers,
                                        auth=self.auth,
                                        verify=self.verify_certificate,
                                        )
            # 502 Bad Gateway triggers retries, since the proxy web
            # server (eg Nginx or Apache) in front of MyTardis could be
            # temporarily restarting
            if response.status_code == 502:
                self._raise_request_exception(response)

        except requests.exceptions.RequestException as e:
            logging.error("Request failed : %s : %s", e.message, url)
            raise e

        return response

    def _md5_python(self, file_path, blocksize=None):
        """
        Calculates the MD5 checksum of a file, returns the hex digest as a
        string. Streams the file in chunks of 'blocksize' to prevent running
        out of memory when working with large files.

        :type file_path: string
        :type blocksize: int
        :return: The hex encoded MD5 checksum.
        :rtype: str
        """
        if not blocksize:
            blocksize = 128

        md5 = hashlib.md5()
        with open(file_path, 'rb') as f:
            while True:
                chunk = f.read(blocksize)
                if not chunk:
                    break
                md5.update(chunk)
        return md5.hexdigest()

    def _md5_subprocess(self, file_path,
                        md5sum_executable='/usr/bin/md5sum'):
        """
        Calculates the MD5 checksum of a file, returns the hex digest as a
        string. Streams the file in chunks of 'blocksize' to prevent running
        out of memory when working with large files.

        :type file_path: string
        :return: The hex encoded MD5 checksum.
        :rtype: str
        """
        out = subprocess.check_output([md5sum_executable, file_path])
        checksum = out.split()[0]
        if len(checksum) == 32:
            return checksum
        else:
            raise ValueError('md5sum failed: %s', out)

    def _md5_file_calc(self, file_path, blocksize=None,
                       subprocess_size_threshold=10*1024*1024,
                       md5sum_executable='/usr/bin/md5sum'):
        """
        Calculates the MD5 checksum of a file, returns the hex digest as a
        string. Streams the file in chunks of 'blocksize' to prevent running
        out of memory when working with large files.

        If the file size is greater than subprocess_size_threshold and the
        md5sum tool exists, spawn a subprocess and use 'md5sum', otherwise
        use the native Python md5 method (~ 3x slower).

        :type file_path: string
        :type blocksize: int
        :param subprocess_size_threshold: Use the md5sum tool via a subprocess
                                           for files larger than this. Otherwise
                                           use Python native method.
        :type subprocess_size_threshold: int
        :return: The hex encoded MD5 checksum.
        :rtype: str
        """
        if os.path.getsize(file_path) > subprocess_size_threshold and \
                os.path.exists(md5sum_executable) and \
                os.access(md5sum_executable, os.X_OK):
            return self._md5_subprocess(file_path,
                                        md5sum_executable=md5sum_executable)
        else:
            return self._md5_python(file_path, blocksize=blocksize)

    def _send_datafile(self, data, filename=None):
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

            response = self.do_post_request('dataset_file',
                                            form,
                                            extra_headers=headers)
            return response

    def _register_datafile_staging(self, data):
        raise NotImplementedError("Registering datafiles in a staging location"
                                  "is not currently implemented.")

    def _register_datafile_shared_storage(self, data):
        response = self.do_post_request('dataset_file', data)
        return response

    def _get_path_from_url(self, url_string):
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

        # logging.debug("Experiment JSON: %s", expt_json)

        data = self.do_post_request('experiment', expt_json)

        if not data.ok or 'Location' not in data.headers:
            logging.error("Creating experiment failed: %s (%s) - %s",
                         data.reason,
                         data.status_code,
                         data.content)
            sys.exit(1)

        return data.headers.get('Location', None)

    def query_instrument(self, name):
        query_params = {u'name': name}
        response = self.do_get_request('instrument', query_params)
        return response.json()

    def query_group(self, name):
        query_params = {u'name': name}
        response = self.do_get_request('group', query_params)
        return response.json()

    def query_user(self, name):
        query_params = {u'username': name}
        response = self.do_get_request('user', query_params)
        return response.json()

    def query_objectacl(self, object_id,
                        content_type='experiment',
                        acl_ownership_type=u'Owner-owned'):

        acl_ownership_type = self._get_ownership_int(acl_ownership_type)

        query_params = {u'object_id': object_id,
                        u'content_type': content_type,
                        u'aclOwnershipType': acl_ownership_type}

        response = self.do_get_request('objectacl',
                                       query_params)
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
            instruments_returned = self.query_instrument(instrument)['objects']
            if instruments_returned:
                instrument_resource_uri = instruments_returned[0]['resource_uri']
                dataset['instrument'] = instrument_resource_uri
            else:
                logging.warning('Querying instrument definition %s failed - is '
                                'the Instrument record defined on the server ?'
                               % instrument)

        dataset_json = MyTardisUploader.dict_to_json(dataset)
        data = self.do_post_request('dataset', dataset_json)

        if not data.ok or 'Location' not in data.headers:
            logging.error("Creating dataset failed: %s", data.text)
            sys.exit(1)

        return data.headers.get('Location', None)

    def should_exclude(self, filepath):
        """
        Returns True if the file path matches one of the patterns to exclude.

        :param filepath:
        :type filepath: str
        :return:
        :rtype: bool
        """
        for p in self.exclude_patterns:
            if p.match(filepath):
                return True
        return False

    def _prepare_file_upload_object(self, file_path, dataset_url_path,
                                    parameter_sets_list=None,
                                    replica_url='',
                                    md5_checksum=None,
                                    base_dir=None):

        if not parameter_sets_list:
            parameter_sets_list = []

        filename = os.path.basename(file_path)
        file_path = os.path.normpath(file_path)

        if base_dir is None:
            base_dir = os.path.dirname(file_path)

        if not replica_url and self.storage_mode == 'shared':
            # replica_url = os.path.normpath(os.path.relpath(file_path,
            #                                start=self.storage_box_location))
            replica_url = file_path.lstrip(self.storage_box_location)
            replica_url = replica_url.lstrip(os.sep)

        replica_list = [{u'url': replica_url,
                         u'location': self.storage_box_name,
                         u'protocol': u'file'},
                        ]
        file_size = os.path.getsize(file_path)
        # Hack to work around MyTardis not accepting
        # files of zero bytes
        # file_size = (file_size if file_size > 0 else -1)
        if md5_checksum is None and not self.fast_mode:
            md5_checksum = self._md5_file_calc(file_path)
        else:
            md5_checksum = '__undetermined__'

        directory = os.path.dirname(file_path).lstrip(os.sep)
        if is_in_tree(file_path, base_dir):
            directory = os.path.normpath(
                os.path.relpath(
                    os.path.dirname(file_path), start=base_dir)).lstrip(os.sep)
        if directory == '.':
            directory = ''

        file_dict = {
            u'dataset': dataset_url_path,
            u'filename': filename,
            u'directory': directory,
            u'md5sum': md5_checksum,
            u'mimetype': mimetypes.guess_type(file_path)[0],
            u'size': file_size,
            u'parameter_sets': parameter_sets_list,
            u'replicas': replica_list,
        }

        return file_dict

    def upload_file(self, file_path, dataset_url_path,
                    parameter_sets_list=None,
                    replica_url='',
                    md5_checksum=None,
                    base_dir=None):

        if self.should_exclude(file_path):
            logging.debug("Skipping %s, matches an exclude pattern.", file_path)
            return None

        file_dict = self._prepare_file_upload_object(
            file_path,
            dataset_url_path,
            parameter_sets_list=parameter_sets_list,
            replica_url=replica_url,
            md5_checksum=md5_checksum,
            base_dir=base_dir)

        if self.storage_mode == 'shared':
            resp_data = self._register_datafile_shared_storage(
                self.dict_to_json(file_dict))
        elif self.storage_mode == 'staging':
            resp_data = self._register_datafile_staging(
                self.dict_to_json(file_dict))
        elif self.storage_mode == 'upload':
            # file_dict.pop(u'replicas', None)
            resp_data = self._send_datafile(
                self.dict_to_json(file_dict),
                filename=file_path)
        else:
            # we should never get here
            raise Exception("Invalid storage mode: " + self.storage_mode)

        if not resp_data.ok or 'Location' not in resp_data.headers:
            logging.error("Registration of data file failed: %s",
                          resp_data.text)
            sys.exit(1)

        return resp_data.headers.get('Location', None)

    def _resource_uri_to_id(self, uri):
        """
        Takes resource URI like: http://example.org/api/v1/experiment/998
        and returns just the id value (998).

        :type uri: str
        :rtype: int
        """
        resource_id = int(urlparse(uri).path.rstrip(os.sep).split(os.sep).pop())
        return resource_id

    def _get_ownership_int(self, ownership_type):
        ownership_type_mappings = {
            u'Owner-owned': 1,
            u'System-owned': 2,
        }
        v = ownership_type_mappings.get(ownership_type, None)
        if v is None:
            raise ValueError("Valid values of acl_ownership_type are %s" %
                             ' or '.join(ownership_type_mappings.keys()))
        return v

    def _share_experiment(self,
                          content_object,
                          plugin_id,
                          entity_id,
                          content_type=u'experiment',
                          acl_ownership_type=u'Owner-owned'):
        """
        Executes an HTTP request to share an MyTardis object with a user or
        group, via updating the ObjectACL.

        :param content_object: The integer ID or URL path to the Experiment,
                               Dataset or DataFile to update.
        :param plugin_id: django_user or django_group
        :param content_type: Django ContentType for the target object, usually
                             'experiment', 'dataset' or 'datafile'
        :type content_object: union(str, int)
        :type plugin_id: str
        :type content_type: basestring
        :return: A requests Response object
        :rtype: Response
        """
        if isinstance(content_object, six.string_types):
            object_id = self._resource_uri_to_id(content_object)
        elif isinstance(content_object, int):
            object_id = content_object
        else:
            raise TypeError("'content_object' must be a URL string or int ID")

        acl_ownership_type = self._get_ownership_int(acl_ownership_type)

        data = {
            u'pluginId': plugin_id,
            u'entityId': str(entity_id),
            u'content_type': str(content_type),
            u'object_id': str(object_id),
            u'aclOwnershipType': acl_ownership_type,
            u'isOwner': True,
            u'canRead': True,
            u'canWrite': True,
            u'canDelete': False,
            u'effectiveDate': None,
            u'expiryDate': None
        }

        response = self.do_post_request('objectacl', self.dict_to_json(data))
        return response

    def share_experiment_with_group(self, experiment, group_name,
                                    *args, **kwargs):
        """
        Executes an HTTP request to share an experiment with a group,
        via updating the ObjectACL.

        :param experiment: The integer ID or URL path to the Experiment.
        :param group_name: Name of the group we with share with.
        :type experiment: union(str, int)
        :type group_name: str
        :return: A requests Response object
        :rtype: Response
        """
        q = self.query_group(group_name)['objects']
        if q:
            group_id = q.pop()['id']
        else:
            raise ValueError("Group not found: %s" % group_name)

        return self._share_experiment(experiment,
                                      'django_group',
                                      group_id,
                                      *args,
                                      **kwargs)

    def share_experiment_with_user(self, experiment, username, *args, **kwargs):
        """
        Executes an HTTP request to share an experiment with a user,
        via updating the ObjectACL.

        :param experiment: The integer ID or URL path to the Experiment
        :param username: The username of the User to share with.
        :type experiment: union(str, int)
        :type username: str
        :return: A requests Response object
        :rtype: Response
        """

        user_id = self.query_user(username)['objects'][0]['id']
        return self._share_experiment(experiment,
                                      'django_user',
                                      user_id,
                                      *args,
                                      **kwargs)

    def remove_object_acl(self, object_id):
        self.do_get_request('objectacl', {})


def setup_logging(logging_config_file='logging_config.toml'):
    import logging

    if logging_config_file and os.path.exists(logging_config_file):
        import logging.config

        # logging.config.fileConfig(logging_config_file)

        import pytoml as toml
        with open(logging_config_file) as f:
            logging_config = toml.load(f)
        logging.config.dictConfig(logging_config)

    else:
        logger = logging.getLogger()
        logger.setLevel(logging.DEBUG)

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

        low_level_request_logger = False
        if low_level_request_logger:
            try:
                import http.client as http_client
            except ImportError:
                # Python 2
                import httplib as http_client

            http_client.HTTPConnection.debuglevel = 1

            requests_log = logging.getLogger("requests.packages.urllib3")
            requests_log.setLevel(logging.DEBUG)
            requests_log.propagate = True

        return logging.getLogger()


def setup_commandline_args(parser=None):
    """
    Takes an argparse.ArgumentParser-like object and adds commandline
    parameters to detect and capture values from.
    :type parser: ArgumentParser
    :return:
    """
    if parser is None:
        parser = ArgumentParser()

    parser.add_argument('--config',
                        dest="config_file",
                        type=str,
                        metavar="MYTARDIS_UPLOADER_CONFIG")
    parser.add_argument("--logging-config",
                        dest="logging_config",
                        type=str,
                        default="logging_config.toml",
                        help="The path to the logging config file "
                             "eg logging_config.toml",
                        metavar="LOGGING_CONFIG")
    parser.add_argument("-f", "--path",
                        dest="path",
                        type=str,
                        help="The absolute PATH of the experiment to be "
                             "uploaded (eg /zfs/storage/awesome_expt_001)",
                        metavar="PATH")
    parser.add_argument("--storage-base-path",
                        dest="storage_base_path",
                        type=str,
                        help="The STORAGE_BASE_PATH of all experiments,"
                             "when using 'shared' storage mode. This should be"
                             "an absolute path, and a prefix of the path"
                             "provided with --path (eg /zfs/storage/)",
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
                             "the command line. Use --api-key instead (and "
                             "ideally only ever from a read-only config file).",
                        metavar="PASSWORD")
    parser.add_argument('--api-key',
                        dest='api_key',
                        type=str,
                        help="A MyTardis (tastypie) REST API key, provided by "
                             "the server.",
                        metavar='API_KEY')
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

    return parser


def parse_settings_yaml(config_filepath):
    with open(config_filepath, 'r') as f:
        parser = SettingsParser(yaml_file=f)

    return parser


def get_config(preparser=None, default_config_filename='uploader_config.yaml'):
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

    :type default_config_filename: str
    :return: The parser.
    :rtype: argparse.ArgumentParser
    """

    if preparser is None:
        preparser = ArgumentParser()

    setup_commandline_args(preparser)
    precheck_options, argv = preparser.parse_known_args()

    parser = None

    if precheck_options.config_file:
        try:
            parser = parse_settings_yaml(precheck_options.config_file)
        except IOError:
            preparser.error("Cannot read config file: %s" %
                            precheck_options.config_file)
    elif os.path.isfile(default_config_filename):
        parser = parse_settings_yaml(default_config_filename)
    else:
        parser = SettingsParser()

    setup_commandline_args(parser)

    return parser


def validate_config(parser, options):
    """
    Validates config options, throws errors and stops excution if
    there is an issue (invalid value or required value missing).

    :rtype : object
    :type options: object
    :type parser: argparse.ArgumentParser
    :return:
    """

    if not options.path:
        parser.error('File path not given (--path)')

    if options.path and not os.path.isabs(options.path):
        parser.error('Path must be an absolute path (--path)')

    if not options.url:
        parser.error('URL to MyTardis instance not given (--url)')

    if not options.username:
        parser.error('MyTardis username not given (--username)')

    valid_storage_modes = ['upload', 'staging', 'shared']
    if options.storage_mode and \
       options.storage_mode not in valid_storage_modes:
        parser.error('--storage-mode must be one of: ' +
                     ', '.join(valid_storage_modes))

    if options.storage_mode == 'shared' and not options.storage_base_path:
        parser.error("--storage-base-path (storage_base_path) must be "
                     "specified when using 'shared' storage mode.")

    if options.storage_mode == 'shared' and \
            not os.path.isabs(options.storage_base_path):
        parser.error("--storage-base-path must be an absolute path "
                     "when using 'shared' storage mode.")

    # if options.storage_mode == 'shared' and \
    #     not is_in_tree(options.storage_base_path, options.path):
    #     parser.error("--path must be in --storage-base-path when using "
    #                  "'shared' storage mode.")

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

    :type exclude_patterns: list[str]
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

    return exclude_regexes


def run():
    ####
    # Le Script
    ####
    # steve.androulakis@monash.edu
    ####
    import getpass

    print("""\
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

    """)

    parser = ArgumentParser()
    parser = setup_commandline_args(parser)
    commandline_options = parser.parse_args()
    config_file = commandline_options.config_file

    try:
        config_options = config_helper.get_config_toml(config_file)
    except IOError as e:
        parser.error("Cannot read config file: %s" % config_file)
        sys.exit(1)

    # Any option that is unset in the commandline args will
    # receive a value from a config file option
    parser.set_defaults(**config_options)

    options = parser.parse_args()

    setup_logging(logging_config_file=options.logging_config)

    validate_config(parser, options)

    password = options.password
    if not password and not options.api_key:
        password = getpass.getpass()

    if password:
        logging.warning("Using 'password' rather than 'api_key' - this is less "
                        "secure and not encouraged.")

    file_path = options.path
    if file_path is '.':
        file_path = os.getcwd()
    title = options.title
    institute = options.institute
    description = options.description
    test_run = options.dry_run

    if options.exclude:
        logging.info("Ignoring files that match: %s\n",
                     ' | '.join(options.exclude))

    mytardis_uploader = MyTardisUploader(
        options.url,
        options.username,
        password=password,
        api_key=options.api_key,
        storage_mode=options.storage_mode,
        storage_box_location=options.storage_base_path,
        storage_box_name=options.storage_box_name,
        verify_certificate=getattr(options, 'verify_certificate', True),
        fast_mode=options.fast,
        exclude_patterns=options.exclude)

    mytardis_uploader.upload_directory(
        file_path,
        title=title,
        description=description,
        institute=institute,
        test_run=test_run)

if __name__ == "__main__":
    run()
