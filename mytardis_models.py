
class MyTardisModelBase(object):
    def __init__(self):
        self.parameter_sets = []

    def fields_to_dict(self, ignore_none=True):
        """
        Return the attributes of the class as a dictionary, ignoring
        the parameter_sets attribute and anything beginning with an underscore.

        These values represent the data for a MyTardis Django model,
        eg and Experiment, Dataset or DataFile, without the additional
        Parameters/ParameterNames.

        :param ignore_none:
        :type ignore_none: bool
        :rtype: dict
        """
        fields = {}
        for k, v in self.__dict__.items():
            if not k.startswith('_') and k != 'parameter_sets':
                if ignore_none:
                    if v is not None:
                        fields[k] = v
                else:
                    fields[k] = v

        return fields

    def package(self):
        request_dict = self.fields_to_dict()
        request_dict['parameter_sets'] = self.package_parameter_sets()
        return request_dict

    @property
    def parameters(self):
        """
        Get just the first parameter set.

        :rtype: MyTardisParameterSet
        """
        if len(self.parameter_sets) >= 1:
            return self.parameter_sets[0]
        else:
            return []

    @parameters.setter
    def parameters(self, value):
        """
        Assign a single parameter set. Deletes any existing parameter sets.

        :type value: MyTardisParameterSet
        """
        del self.parameter_sets[:]
        self.parameter_sets = [value]

    @parameters.deleter
    def parameters(self):
        """
        Deletes all parameter sets.
        """
        del self.parameter_sets[:]

    def package_parameter_sets(self):
        """

        Returns:
        [{'schema':'http://x', 'parameters': [{'name':'w', 'value': 'z'},
         {'schema':'http://y', 'parameters': [{'name':'a', 'value': 'b'},
        ]

        :return:
        :rtype: list[dict]
        """
        parameter_sets = [param_set.package_parameter_set()
                          for param_set in self.parameter_sets]

        return parameter_sets


# TODO: Make this the parent class of parameter sets
# (eg IlluminaSequencingRunBase) and add a parameter_sets = []
# as an attribute to MyTardisModel. This way we can easily keep a
# separation between parameter sets and attributes on the model itself
# (eg 'title'). MyTardisModel (and subclasses such as Experiment) should
# be able to generate the full dictionary required for an API call,
# eg {'title': 'bla', 'description': 'foo',
#       'parameter_sets': [{'schema':'http://x',
#                           'parameters': [{'name':'x', 'value': 'y'},]
#
class MyTardisParameterSet(object):
    def __init__(self):
        self._namespace__schema = None  # type: str

    def from_dict(self, d, existing_only=True):
        for k, v in d.items():
            if existing_only:
                if hasattr(self, k):
                    setattr(self, k, v)
                else:
                    raise KeyError('Attribute %s not found in %s' %
                                   (k, self.__class__))
            else:
                setattr(self, k, v)

    def to_dict(self, ignore_none=True):
        params = {}
        for k, v in self.__dict__.items():
            if not k.startswith('_'):
                if ignore_none:
                    if v is not None:
                        params[k] = v
                else:
                    params[k] = v

        return params

    def package_parameter_list(self, ignore_none=True):
        pset = self.to_dict(ignore_none=ignore_none)
        parameter_list = [{u'name': k, u'value': v} for k, v in pset.items()]
        return parameter_list

    def package_parameter_set(self, ignore_none=True):
        return {u'schema': self._namespace__schema,
                u'parameters': self.package_parameter_list(
                    ignore_none=ignore_none)
                }

    def to_schema(self):
        schema = {'pk': None, 'model': '', 'fields': {}}
        for k, v in self.__dict__.items():
            if k.startswith('_') and k.endswith('__schema'):
                if k in ['__pk_schema', '__model_schema']:
                    schema[k[2:]] = v
                else:
                    schema['fields'][k[1:]] = v

        return [schema]

    def to_parameter_schema(self):
        param_schemas = []
        for k, v in self.__dict__.items():
            if k.startswith('_') and k.endswith('__attr_schema'):
                param_schemas.append(v)

        return param_schemas


class Experiment(MyTardisModelBase):
    def __init__(self, *args, **kwargs):
        super(Experiment, self).__init__(*args, **kwargs)
        self.title = None  # type: str
        self.institution_name = None  # type: str
        self.description = None  # type: str
        self.start_time = None  # type: datetime
        self.end_time = None  # type: datetime
        self.created_time = None  # type: datetime
        self.created_by = None  # type: str
        self.handle = None  # type: str
        self.locked = None  # type: bool
        self.public_access = None  # type: int
        self.license = None  # type: str


class Dataset(MyTardisModelBase):
    def __init__(self, *args, **kwargs):
        super(Dataset, self).__init__(*args, **kwargs)
        pass


class DataFile(MyTardisModelBase):
    def __init__(self, *args, **kwargs):
        super(DataFile, self).__init__(*args, **kwargs)
        pass
