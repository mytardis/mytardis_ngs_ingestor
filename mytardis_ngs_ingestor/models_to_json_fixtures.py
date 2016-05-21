#!/usr/bin/env python
# Serializes MyTardisParameterSet models to JSON fixtures, ready to populate
# the MyTardis database (eg via ./mytardis.py loaddata)

# Usage:
# $ models_to_json_fixtures.py >sequencing_facility_fixtures.json

from __future__ import print_function

import illumina_uploader

if __name__ == "__main__":
    print(illumina_uploader.dump_schema_fixtures_as_json())

