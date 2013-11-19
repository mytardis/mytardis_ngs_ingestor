# MyTardis Uploader v1.1
# Steve Androulakis <steve.androulakis@monash.edu>
# Thanks Grischa Meyer <grischa.meyer@monash.edu> for initial script

import urllib2
import base64
import os
import mimetypes
import json
import requests
from requests.auth import HTTPBasicAuth
from time import strftime


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
                        test_run=False
                        ):

        title = title or os.path.basename(os.path.abspath(file_path))
        print 'Creating experiment: %s' % title

        created = False
        exp_url = "/test/"
        if not test_run:
            exp_url = self.create_experiment(title,\
                                    institute or '',
                                    '%s \n\nAutomatically generated on %s by https://github.com/steveandroulakis/mytardis-uploader'\
                                       % (description or 'No description.', strftime("%Y-%m-%d %H:%M:%S")))
            created = True

        for item in os.listdir(file_path):

            if item.startswith('.'):
                continue  # filter files/dirs starting with .

            if os.path.isfile(os.path.join(file_path, item)):
                print 'Skipping root-level file: %s' % item
            else:
                print '\tFound directory: %s' % item

                print '\tCreating dataset: %s' % item

                ds_url = "/test/"
                if not test_run:
                    ds_url = self.create_dataset('%s' % item,
                              [self._get_path_from_url(exp_url)],
                              )

                for dirname, dirnames, filenames in os.walk(os.path.join(file_path, item)):

                    for filename in filenames:

                        if filename.startswith('.'):
                            continue  # filter files/dirs starting with .

                        sub_file_path = '%s/%s' % (dirname, filename)
                        print "\t\tUploading file '%s' to dataset '%s'." % (sub_file_path, item)

                        #f_url = "/test/"
                        if not test_run:
                            self.upload_file(sub_file_path, self._get_path_from_url(ds_url))
                            #print f_url

        if created:
            exp_id = exp_url.rsplit('/')[-2]
            new_exp_url = "%s/experiment/view/%s/" % (self.mytardis_url, exp_id)
            print "Experiment created: %s" % new_exp_url
            return new_exp_url

        else:
            print "Dry run complete."
            return "http://example.com/test/success"

    def _send_data(self, data, urlend, method="POST"):
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
        output = "error"
        output = urllib2.urlopen(myrequest)

        return output

    def _md5_file_calc(self, file_path):
        import hashlib
        return hashlib.md5(open(file_path, 'rb').read()).hexdigest()

    def _send_datafile(self, data, urlend, method='POST', filename=None):
        url = self.v1_api_url % urlend
        headers = {'Accept': 'application/json'}
        #import ipdb; ipdb.set_trace()
        response = requests.post(url, data={"json_data": data}, headers=headers,
                                 files={'attached_file': open(filename, 'rb')},
                                 auth=HTTPBasicAuth(self.username, self.password)
                                 )
        # for item in response:
        #     print item
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
    # parameter_dict.append({u'name': 'diffractometerType', u'value': 'this is my value'})
        parameter_set = {}
        parameter_set['schema'] = schema
        parameter_set['parameters'] = parameter_list
        return parameter_set

    def create_experiment(self, title, institution, description):

        exp_dict = {
            u'description': description,
            u'institution_name': institution,
            u'title': title
            }

        exp_json = json.dumps(exp_dict)

        data = self._send_data(exp_json, 'experiment/')

        return data.info().getheaders('Location')[0]

    def create_dataset(self, description, experiments_list, parameter_sets_list=[], immutable=False):

        dataset_dict = {
                          u'description': description,
                          u'experiments': experiments_list,
                          u'immutable': immutable,
                          u'parameter_sets': parameter_sets_list
                        }

        dataset_json = json.dumps(dataset_dict)

        data = self._send_data(dataset_json, 'dataset/')

        return data.info().getheaders('Location')[0]

    def upload_file(self, file_path, dataset_path):
    #print upload_file('cli.py', '/api/v1/dataset/143/').headers['location']

        file_dict = {
            u'dataset': dataset_path,
            u'filename': os.path.basename(file_path),
            u'md5sum': self._md5_file_calc(file_path),
            u'mimetype': mimetypes.guess_type(file_path)[0],
            u'size': os.path.getsize(file_path)
            }

        file_json = json.dumps(file_dict)
        data = self._send_datafile(file_json, 'dataset_file/', 'POST', file_path)

        return data.headers['location']


def run():
    ####
    # Le Script
    ####
    #   steve.androulakis@monash.edu
    ####

    from optparse import OptionParser
    import getpass

    print "MyTardis uploader generic v1"
    print "Steve Androulakis <steve.androulakis@monash.edu>"
    print "Uploads the given directory as an Experiment, and the immediate \n" \
          "sub-directories below it as Datasets in MyTardis."
    print "Eg. python mytardis_uploader.py -l http://mytardis-server.com.au -u steve" \
          " -f /Users/steve/Experiment1/"
    print ""

    parser = OptionParser()
    parser.add_option("-f", "--path", dest="file_path",
                      help="The PATH of the experiment to be uploaded", metavar="PATH")
    parser.add_option("-l", "--url", dest="mytardis_url",
                      help="The URL to the MyTardis installation", metavar="URL")
    parser.add_option("-u", "--username", dest="username",
                      help="Your MyTardis USERNAME", metavar="USERNAME")
    parser.add_option("-t", "--title", dest="title",
                      help="Experiment TITLE", metavar="TITLE")
    parser.add_option("-d", "--description", dest="description",
                      help="Experiment DESCRIPTION", metavar="DESCRIPTION")
    parser.add_option("-i", "--institute", dest="institute",
                      help="Experiment INSTITUTE (eg university)", metavar="INSTITUTE")
    parser.add_option("-r", "--dry",
                      action="store_true", dest="dry_run", default=False,
                      help="Dry run (don't create anything)")

    (options, args) = parser.parse_args()

    if not options.file_path:
        parser.error('file path not given')

    if not options.mytardis_url:
        parser.error('url to MyTardis not given')

    if not options.username:
        parser.error('MyTardis username not given')

    pw = getpass.getpass()

    file_path = options.file_path
    title = options.title
    institute = options.institute
    description = options.description
    test_run = options.dry_run
    mytardis_url = options.mytardis_url
    username = options.username
    password = pw

    mytardis_uploader = MyTardisUploader(mytardis_url,
                                         username,
                                         password)

    mytardis_uploader.upload_directory(file_path,
                                       title=title,
                                       description=description,
                                       institute=institute,
                                       test_run=test_run)

if __name__ == "__main__":
    run()
