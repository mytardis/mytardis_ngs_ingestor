# Based on: https://fadeit.dk/blog/2015/10/12/mailgun-python-log-handler/
import os
import logging
import requests


class MailgunHandler(logging.Handler):

    def __init__(self, api_url, api_key, sender, recipients, subject_prefix):
        # run the regular Handler __init__
        logging.Handler.__init__(self)
        self.api_url = api_url
        self.api_key = api_key
        self.sender = sender
        self.recipients = recipients
        self.hostname = os.uname()[1]
        self.subject_prefix = subject_prefix

    def _format_subject(self, message):
        return "%s (%s): %s" % (self.subject_prefix,
                                self.hostname,
                                message)

    def emit(self, record):
        # record.message is the log message (only created after self.format is
        # called)
        for recipient in self.recipients:
            text = self.format(record)
            subject = self._format_subject(record.message)
            data = {
                "from": self.sender,
                "to": recipient,
                "subject": subject,
                "text": text,
            }
            requests.post(self.api_url, auth=("api", self.api_key), data=data)
