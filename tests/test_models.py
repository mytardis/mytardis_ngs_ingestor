import os
import unittest
from mytardis_ngs_ingestor.mytardis_models import MyTardisParameterSet


class TestingParameterSetModel(MyTardisParameterSet):
    """
    :type some_string_field: unicode
    :type a_float_field: float
    """

    def __init__(self):
        super(TestingParameterSetModel, self).__init__()
        self.some_string_field = None  # type: unicode
        self.a_float_field = None  # type: float

        # Dictionaries to allow reconstitution of the schema for each parameter

        # run_id fixture
        self._run_id__attr_schema = {
            u'pk': None,
            u'model': u'tardis_portal.parametername',
            u'fields': {u'name': u'some_string_field',
                        u'data_type': 2,
                        u'immutable': True,
                        u'is_searchable': True,
                        u'choices': u'',
                        u'comparison_type': 1,
                        u'full_name': u'A string field',
                        u'units': u'', u'order': 9999,
                        u'schema': [u'http://www.tardis.edu.au/schema/test']}}  # type: dict

        # run_number fixture
        self._run_number__attr_schema = {
            u'pk': None,
            u'model': u'tardis_portal.parametername',
            u'fields': {u'name': u'a_float_field',
                        u'data_type': 1,
                        u'immutable': True,
                        u'is_searchable': True,
                        u'choices': u'',
                        u'comparison_type': 1,
                        u'full_name': u'A float field',
                        u'units': u'',
                        u'order': 9999,
                        u'schema': [u'http://www.tardis.edu.au/schema/test']}}  # type: dict

        self._subtype__schema = "testing-parameter-set"  # type: unicode
        self._model__schema = "tardis_portal.schema"  # type: unicode
        self._name__schema = "A test parameter set model"  # type: unicode
        self._pk__schema = None  # type: NoneType
        self._type__schema = 1  # type: int
        self._hidden__schema = False  # type: bool
        self._namespace__schema = "http://www.tardis.edu.au/schema/test"  # type: unicode
        self._immutable__schema = True  # type: bool


class IlluminaParserTestCase(unittest.TestCase):
    def setUp(self):
        self.parameterset_model = TestingParameterSetModel()

    def tearDown(self):
        pass

    def test__to_schema(self):
        schema_as_dict = self.parameterset_model.to_schema()
        self.assertEqual(schema_as_dict, [
            {'pk': None, 'model': 'tardis_portal.schema',
             'fields': {'name': 'A test parameter set model',
                        'namespace': 'http://www.tardis.edu.au/schema/test',
                        'subtype': 'testing-parameter-set', 'hidden': False,
                        'type': 1, 'immutable': True}}])

    def test__to_parameter_schema(self):
        param_schema_dict = self.parameterset_model.to_parameter_schema()
        self.maxDiff = None
        self.assertDictEqual(param_schema_dict[0],
                          {u'pk': None,
                           u'model': u'tardis_portal.parametername',
                           u'fields': {u'full_name': u'A float field',
                                       u'comparison_type': 1,
                                       u'schema': [
                                       u'http://www.tardis.edu.au/schema/test'],
                                       u'name': u'a_float_field',
                                       u'data_type': 1,
                                       u'units': u'', u'order': 9999,
                                       u'immutable': True,
                                       u'is_searchable': True,
                                       u'choices': u''}})
        self.assertDictEqual(param_schema_dict[1],
                          {u'pk': None,
                           u'model': u'tardis_portal.parametername',
                           u'fields': {u'full_name': u'A string field',
                                       u'comparison_type': 1,
                                       u'schema': [
                                       u'http://www.tardis.edu.au/schema/test'],
                                       u'name': u'some_string_field',
                                       u'data_type': 2,
                                       u'units': u'', u'order': 9999,
                                       u'immutable': True,
                                       u'is_searchable': True,
                                       u'choices': u''}})


if __name__ == '__main__':
    unittest.main()
