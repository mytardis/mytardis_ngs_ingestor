import re


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

    def package(self, ignore_none=True):
        request_dict = self.fields_to_dict(ignore_none=ignore_none)
        request_dict['parameter_sets'] = \
            self.package_parameter_sets(ignore_none=ignore_none)
        return request_dict

    def to_json(self, ignore_none=True):
        import datetime, json
        date_handler = lambda obj: (
            obj.isoformat()
            if isinstance(obj, datetime.datetime)
            or isinstance(obj, datetime.date)
            else None
        )
        return json.dumps(self.package(ignore_none=ignore_none),
                          default=date_handler)

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

    def package_parameter_sets(self, ignore_none=True):
        """

        Returns:
        [{'schema':'http://x', 'parameters': [{'name':'w', 'value': 'z'},
         {'schema':'http://y', 'parameters': [{'name':'a', 'value': 'b'},
        ]

        :return:
        :rtype: list[dict]
        """
        parameter_sets = [param_set.package_parameter_set(ignore_none=
                                                          ignore_none)
                          for param_set in self.parameter_sets]

        return parameter_sets


class MyTardisParameterSet(object):
    def __init__(self):
        self._namespace__schema = None  # type: unicode

    def from_dict(self, d, existing_only=True, ignore_missing=True):
        for k, v in d.items():
            if existing_only:
                if hasattr(self, k):
                    setattr(self, k, v)
                else:
                    if not ignore_missing:
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
        schema = {'pk': None, 'model': 'tardis_portal.schema', 'fields': {}}
        attributes = sorted(self.__dict__.items())
        for k, v in attributes:
            if k.startswith('_') and k.endswith('__schema'):
                # remove leading _ and __schema suffix
                kname = re.sub(re.escape('__schema')+'$', '', k[1:])
                # model and pk go into the top level, all the rest are
                # part of the nest 'fields' dictionary
                if kname in ['pk', 'model']:
                    schema[kname] = v
                else:
                    schema['fields'][kname] = v

        return [schema]

    def to_parameter_schema(self):
        param_schemas = []
        attributes = sorted(self.__dict__.items())
        for k, v in attributes:
            if k.startswith('_') and k.endswith('__attr_schema'):
                param_schemas.append(v)

        param_schemas.sort()
        return param_schemas


class Experiment(MyTardisModelBase):
    def __init__(self, *args, **kwargs):
        super(Experiment, self).__init__(*args, **kwargs)
        self.title = None  # type: unicode
        self.institution_name = None  # type: unicode
        self.description = None  # type: unicode
        self.start_time = None  # type: datetime
        self.end_time = None  # type: datetime
        self.created_time = None  # type: datetime
        self.created_by = None  # type: unicode
        self.handle = None  # type: unicode
        self.locked = None  # type: bool
        self.public_access = None  # type: int
        self.license = None  # type: unicode


class Dataset(MyTardisModelBase):
    def __init__(self, *args, **kwargs):
        super(Dataset, self).__init__(*args, **kwargs)
        self.experiments = []  # type: list[str]
        self.description = None  # type: unicode
        self.immutable = None  # type: bool
        self.instrument = None  # type: unicode


class DataFile(MyTardisModelBase):
    def __init__(self, *args, **kwargs):
        super(DataFile, self).__init__(*args, **kwargs)
        pass
