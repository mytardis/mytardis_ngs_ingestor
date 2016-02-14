#!/usr/bin/python
# Takes a MyTardis fixture defining Schemas and ParameterNames
# (as JSON) and generates a set of Python base classes
# as a data model for use by REST API clients.
#
# Extra 'fields' from the schema are added to classes as underscore prefixed
# _attributes, like '_immutable__schema' and '_immutable__attr_schema'
# so that in theory they can be used to regenerate the original JSON fixtures
# from which they were derived. The goal is to use the class definition as the
# 'one source of truth' and regenerate JSON fixtures when the class is modified.
# Class names are derived from the 'subtype' defined in the schema, and a
# prefixed with 'Base' (eg some-subtype becomes SomeSubtypeBase) - these
# parent classes can be imported and overridden to keep a separation between
# generated code and hand written methods.

#
# Usage:
# python fixture_to_class.py myfixture.json >models.py
#

import sys
import json
import re
from collections import OrderedDict
from textwrap import wrap

INDENT = 4

SCHEMA_TYPES = ('Experiment', 'Dataset', 'DataFile', 'None', 'Instrument')

PARAM_TYPES = (
    ('float', 'NUMERIC'),
    ('unicode', 'STRING'),
    ('unicode', 'URL'),
    ('unicode', 'LINK'),
    ('unicode', 'FILENAME'),
    ('datetime.datetime', 'DATETIME'),
    ('unicode', 'LONGSTRING'),
    ('dict', 'JSON'),
)


# You make wish to modify this to better transform your schema names
# into class names
def sanitize_class_name(name):
    name = name.replace('-', ' ')
    name = name.title()  # to some kind of Camel Case
    # then remove non-alphanumerics and spaces
    name = re.sub('[^0-9a-zA-Z]+', '', name)
    return name


def sanitize_attr_name(name):
    # python attribute names should be lowercase with underscores
    name = name.lower()
    name = re.sub('[^0-9a-zA-Z]+', '_', name)
    return name


def wrap_python_code(code):
    lines = wrap(code, drop_whitespace=False)
    quoted = []
    open_quote_next_line = False
    for l in lines:

        # This code attempts to add missing quotes to strings
        # that are split across lines by 'wrap'. There could
        # be bugs.
        ll = unicode(l)
        unpaired_ticks = ll.count(u"'") % 2
        if open_quote_next_line and unpaired_ticks > 0:
            ll = u"u'" + ll
            open_quote_next_line = False

        unpaired_ticks = ll.count(u"'") % 2
        if unpaired_ticks > 0:
            if ll[-1:] != u"'":
                ll += u"'"
                open_quote_next_line = True
            quoted.append(ll)
        else:
            quoted.append(ll)

    return '\n'.join(quoted)


class ClassDef:
    def __init__(self, name, namespace, parent_class, fixture=None):
        self.name = name
        self.namespace = namespace
        self.parent_class = parent_class
        self.fixture = fixture
        self.attributes = []

    def format_class_attributes(self):
        attribs = []
        for a in self.attributes:
            attribs.append("# %s\n%s = None  # type: %s" %
                           (a.comment,
                            a.name,
                            a.data_type))

        return '\n'.join(attribs)

    def format_instance_attributes(self):
        attribs = []
        for a in self.attributes:
            attribs.append("# %s\nself.%s = None  # type: %s\n" %
                           (a.comment,
                            a.name,
                            a.data_type))

        return '\n'.join(attribs)

    def format_instance_attribute_schema_dicts(self):
        attribs = []
        for a in self.attributes:
            attribs.append(
                "# %s fixture\n" % a.name +
                wrap_python_code("self._%s__attr_schema = %s" %
                                 (a.name, a.fixture)) +
                "  # type: dict\n")

        return '\n'.join(attribs)

    def format_attrib_docstring(self):
        docstr_types = []
        for a in self.attributes:
            docstr_types.append(":type %s: %s" % (a.name, a.data_type))

        return '\n'.join(docstr_types)

    def format_fixture_to_instance_attrs(self):
        attribs = []
        fields = dict(self.fixture['fields'])
        fields['pk'] = self.fixture['pk']
        fields['model'] = self.fixture['model']
        for k, v in fields.items():
            vv = unicode(v)
            if isinstance(v, (str, unicode)):
                vv = '"%s"' % v
            attribs.append("self._%s__schema = %s  # type: %s" %
                           (k, vv, type(v).__name__))

        return '\n'.join(attribs)


class AttributeDef:
    name = None
    data_type = None
    comment = ""
    fixture = None

    def __init__(self, name, data_type, comment='', fixture=None):
        self.name = name
        self.data_type = data_type
        self.comment = comment
        self.fixture = fixture


# http://stackoverflow.com/a/8348914/77990
def indent(lines, amount, ch=' '):
    padding = amount * ch
    return padding + ('\n' + padding).join(lines.split('\n'))


if __name__ == "__main__":
    fixtures = None
    with open(sys.argv[1], 'r') as f:
        fixtures = json.loads(f.read())

    classes = OrderedDict()
    for obj in fixtures:
        fields = obj['fields']
        if obj['model'] == 'tardis_portal.schema':
            name = sanitize_class_name(fields['subtype']) + 'Base'
            schema_type = fields['type'] - 1
            namespace = fields['namespace']
            classes[namespace] = ClassDef(name, namespace,
                                          SCHEMA_TYPES[schema_type],
                                          fixture=obj)
        if obj['model'] == 'tardis_portal.parametername':
            name = sanitize_attr_name(obj['fields']['name'])
            namespace = fields['schema'][0]
            data_type = fields['data_type'] - 1
            type_str = PARAM_TYPES[data_type][0]

            for a in classes[namespace].attributes:
                if name == a.name:
                    sys.stderr.write("WARNING: Duplicate ParameterName: "
                                     "%s (%s)\n" % (name, namespace))

            classes[namespace].attributes.append(
                AttributeDef(name,
                             type_str,
                             comment=fields['full_name'],
                             fixture=obj)
            )

    print "# Data model generated from %s\n\n" % (sys.argv[1])
    print "from mytardis_models import *\n\n"

    for namespace, klass in classes.items():
        print '''
class %(name)s(%(parent)s):
    """
%(docstring_text)s

%(docstring_types)s
    """

    def __init__(self):
        super(%(name)s, self).__init__()
%(self_attributes)s

        # Dictionaries to allow reconstitution of the schema for each parameter

%(parameter_schema_meta)s

%(schema_meta)s
''' % {
            'name': klass.name,
            'parent': 'MyTardisParameterSet',  # klass.parent_class,
            'docstring_text': "",
            # unused
            'fixture_json': indent(json.dumps(klass.fixture, indent=True),
                                   INDENT),
            'docstring_types': indent(klass.format_attrib_docstring(), INDENT),
            'self_attributes': indent(klass.format_instance_attributes(),
                                      INDENT * 2),
            'parameter_schema_meta': indent(
                klass.format_instance_attribute_schema_dicts(), INDENT * 2),
            'schema_meta': indent(klass.format_fixture_to_instance_attrs(),
                                  INDENT * 2)
        }
