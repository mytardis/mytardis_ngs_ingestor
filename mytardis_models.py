
class MyTardisModel(object):
    def __init__(self):
        pass

    def to_schema(self):
        schema = {'pk': None, 'model': '', 'fields': {}}
        for k, v in self.__dict__.items():
            if k[0] == '_' and k[-8:] != '__schema':
                if k in ['_pk', '_model']:
                    schema[k[1:]] = v
                else:
                    schema['fields'][k[1:]] = v

        return [schema]

    def to_parameter_schema(self):
        param_schemas = []
        for k, v in self.__dict__.items():
            if k[0] == '_' and k[-8:] == '__schema':
                    param_schemas.append(v)

        return param_schemas


class Experiment(MyTardisModel):
    def __init__(self):
        super(MyTardisModel, self).__init__()
        pass


class Dataset(MyTardisModel):
    def __init__(self):
        super(MyTardisModel, self).__init__()
        pass


class DataFile(MyTardisModel):
    def __init__(self):
        super(MyTardisModel, self).__init__()
        pass
