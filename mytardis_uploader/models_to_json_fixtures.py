#!/usr/bin/env python
# Serializes MyTardisParameterSet models to JSON fixtures, ready to populate
# the MyTardis database (eg via ./mytardis.py loaddata)

# Usage:
# $ models_to_json_fixtures.py >sequencing_facility_fixtures.json

import json
import models
from mytardis_models import MyTardisParameterSet

if __name__ == "__main__":
    fixtures = []
    for name, klass in models.__dict__.items():
        if (not name.startswith('__') and
                issubclass(klass, MyTardisParameterSet)):
            fixtures += klass().to_schema()
            fixtures += klass().to_parameter_schema()

    print json.dumps(fixtures, indent=2)
