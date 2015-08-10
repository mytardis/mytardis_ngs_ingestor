# Data model generated from ../fixtures/sequencing_facility_schema.json


from mytardis_models import *



class IlluminaSequencingRunBase(Experiment):
    """


    :type run_id: unicode
    :type run_number: float
    :type flowcell_id: unicode
    :type instrument_id: unicode
    :type instrument_model: unicode
    :type read_cycles: unicode
    :type chemistry: unicode
    :type operator_name: unicode
    :type rta_version: unicode
    """

    def __init__(self):
        super(Experiment, self).__init__()
        # Run Unique ID
        self.run_id = None  # type: unicode
        
        # Run number
        self.run_number = None  # type: float
        
        # Flowcell ID
        self.flowcell_id = None  # type: unicode
        
        # Instrument ID
        self.instrument_id = None  # type: unicode
        
        # Instrument model
        self.instrument_model = None  # type: unicode
        
        # Number of cycles in each read [index reads in (brackets)]
        self.read_cycles = None  # type: unicode
        
        # Terminator chemistry
        self.chemistry = None  # type: unicode
        
        # Instrument operator
        self.operator_name = None  # type: unicode
        
        # Illumina RTA version
        self.rta_version = None  # type: unicode
        

        # Dictionaries to allow reconstitution of the schema for each parameter

        # run_id fixture
        self._run_id__schema = {u'pk': None, u'model':
        u'tardis_portal.parametername', u'fields': {u'name': u'run_id',
        u'data_type': 2, u'immutable': True, u'is_searchable': True,
        u'choices': u'', u'comparison_type': 1, u'full_name': u'Run Unique'
        u'ID', u'units': u'', u'order': 9999, u'schema':
        [u'http://www.tardis.edu.au/schemas/ngs/run/illumina']}}  # type: dict
        
        # run_number fixture
        self._run_number__schema = {u'pk': None, u'model':
        u'tardis_portal.parametername', u'fields': {u'name': u'run_number',
        u'data_type': 1, u'immutable': True, u'is_searchable': True,
        u'choices': u'', u'comparison_type': 1, u'full_name': u'Run number',
        u'units': u'', u'order': 9999, u'schema':
        [u'http://www.tardis.edu.au/schemas/ngs/run/illumina']}}  # type: dict
        
        # flowcell_id fixture
        self._flowcell_id__schema = {u'pk': None, u'model':
        u'tardis_portal.parametername', u'fields': {u'name': u'flowcell_id',
        u'data_type': 2, u'immutable': True, u'is_searchable': True,
        u'choices': u'', u'comparison_type': 1, u'full_name': u'Flowcell ID',
        u'units': u'', u'order': 9999, u'schema':
        [u'http://www.tardis.edu.au/schemas/ngs/run/illumina']}}  # type: dict
        
        # instrument_id fixture
        self._instrument_id__schema = {u'pk': None, u'model':
        u'tardis_portal.parametername', u'fields': {u'name': u'instrument_id',
        u'data_type': 2, u'immutable': True, u'is_searchable': True,
        u'choices': u'', u'comparison_type': 1, u'full_name': u'Instrument'
        u'ID', u'units': u'', u'order': 9999, u'schema':
        [u'http://www.tardis.edu.au/schemas/ngs/run/illumina']}}  # type: dict
        
        # instrument_model fixture
        self._instrument_model__schema = {u'pk': None, u'model':
        u'tardis_portal.parametername', u'fields': {u'name':
        u'instrument_model', u'data_type': 2, u'immutable': True,
        u'is_searchable': True, u'choices': u'', u'comparison_type': 1,
        u'full_name': u'Instrument model', u'units': u'', u'order': 9999,
        u'schema': [u'http://www.tardis.edu.au/schemas/ngs/run/illumina']}}  # type: dict
        
        # read_cycles fixture
        self._read_cycles__schema = {u'pk': None, u'model':
        u'tardis_portal.parametername', u'fields': {u'name': u'read_cycles',
        u'data_type': 2, u'immutable': True, u'is_searchable': False,
        u'choices': u'', u'comparison_type': 1, u'full_name': u'Number of'
        u'cycles in each read [index reads in (brackets)]', u'units': u'',
        u'order': 9999, u'schema':
        [u'http://www.tardis.edu.au/schemas/ngs/run/illumina']}}  # type: dict
        
        # chemistry fixture
        self._chemistry__schema = {u'pk': None, u'model':
        u'tardis_portal.parametername', u'fields': {u'name': u'chemistry',
        u'data_type': 2, u'immutable': True, u'is_searchable': False,
        u'choices': u'', u'comparison_type': 1, u'full_name': u'Terminator'
        u'chemistry', u'units': u'', u'order': 9999, u'schema':
        [u'http://www.tardis.edu.au/schemas/ngs/run/illumina']}}  # type: dict
        
        # operator_name fixture
        self._operator_name__schema = {u'pk': None, u'model':
        u'tardis_portal.parametername', u'fields': {u'name': u'operator_name',
        u'data_type': 2, u'immutable': True, u'is_searchable': False,
        u'choices': u'', u'comparison_type': 1, u'full_name': u'Instrument'
        u'operator', u'units': u'', u'order': 9999, u'schema':
        [u'http://www.tardis.edu.au/schemas/ngs/run/illumina']}}  # type: dict
        
        # rta_version fixture
        self._rta_version__schema = {u'pk': None, u'model':
        u'tardis_portal.parametername', u'fields': {u'name': u'rta_version',
        u'data_type': 2, u'immutable': True, u'is_searchable': False,
        u'choices': u'', u'comparison_type': 1, u'full_name': u'Illumina RTA'
        u'version', u'units': u'', u'order': 9999, u'schema':
        [u'http://www.tardis.edu.au/schemas/ngs/run/illumina']}}  # type: dict
        

        self._subtype = "illumina-sequencing-run"  # type: unicode
        self._model = "tardis_portal.schema"  # type: unicode
        self._name = "Illumina Sequencing Run"  # type: unicode
        self._pk = None  # type: NoneType
        self._type = 1  # type: int
        self._hidden = False  # type: bool
        self._namespace = "http://www.tardis.edu.au/schemas/ngs/run/illumina"  # type: unicode
        self._immutable = True  # type: bool


class DemultiplexedSamplesBase(Experiment):
    """


    :type run_id: unicode
    :type run_experiment: unicode
    :type run_number: float
    :type flowcell_id: unicode
    :type instrument_id: unicode
    :type instrument_model: unicode
    :type read_cycles: unicode
    :type chemistry: unicode
    :type operator_name: unicode
    :type rta_version: unicode
    """

    def __init__(self):
        super(Experiment, self).__init__()
        # Run Unique ID
        self.run_id = None  # type: unicode
        
        # Run Experiment link
        self.run_experiment = None  # type: unicode
        
        # Run number
        self.run_number = None  # type: float
        
        # Flowcell ID
        self.flowcell_id = None  # type: unicode
        
        # Instrument ID
        self.instrument_id = None  # type: unicode
        
        # Instrument model
        self.instrument_model = None  # type: unicode
        
        # Number of cycles in each read [index reads in (brackets)]
        self.read_cycles = None  # type: unicode
        
        # Terminator chemistry
        self.chemistry = None  # type: unicode
        
        # Instrument operator
        self.operator_name = None  # type: unicode
        
        # Illumina RTA version
        self.rta_version = None  # type: unicode
        

        # Dictionaries to allow reconstitution of the schema for each parameter

        # run_id fixture
        self._run_id__schema = {u'pk': None, u'model':
        u'tardis_portal.parametername', u'fields': {u'name': u'run_id',
        u'data_type': 2, u'immutable': True, u'is_searchable': True,
        u'choices': u'', u'comparison_type': 1, u'full_name': u'Run Unique'
        u'ID', u'units': u'', u'order': 9999, u'schema':
        [u'http://www.tardis.edu.au/schemas/ngs/project']}}  # type: dict
        
        # run_experiment fixture
        self._run_experiment__schema = {u'pk': None, u'model':
        u'tardis_portal.parametername', u'fields': {u'name':
        u'run_experiment', u'data_type': 4, u'immutable': True,
        u'is_searchable': True, u'choices': u'', u'comparison_type': 1,
        u'full_name': u'Run Experiment link', u'units': u'', u'order': 9999,
        u'schema': [u'http://www.tardis.edu.au/schemas/ngs/project']}}  # type: dict
        
        # run_number fixture
        self._run_number__schema = {u'pk': None, u'model':
        u'tardis_portal.parametername', u'fields': {u'name': u'run_number',
        u'data_type': 1, u'immutable': True, u'is_searchable': True,
        u'choices': u'', u'comparison_type': 1, u'full_name': u'Run number',
        u'units': u'', u'order': 9999, u'schema':
        [u'http://www.tardis.edu.au/schemas/ngs/project']}}  # type: dict
        
        # flowcell_id fixture
        self._flowcell_id__schema = {u'pk': None, u'model':
        u'tardis_portal.parametername', u'fields': {u'name': u'flowcell_id',
        u'data_type': 2, u'immutable': True, u'is_searchable': True,
        u'choices': u'', u'comparison_type': 1, u'full_name': u'Flowcell ID',
        u'units': u'', u'order': 9999, u'schema':
        [u'http://www.tardis.edu.au/schemas/ngs/project']}}  # type: dict
        
        # instrument_id fixture
        self._instrument_id__schema = {u'pk': None, u'model':
        u'tardis_portal.parametername', u'fields': {u'name': u'instrument_id',
        u'data_type': 2, u'immutable': True, u'is_searchable': True,
        u'choices': u'', u'comparison_type': 1, u'full_name': u'Instrument'
        u'ID', u'units': u'', u'order': 9999, u'schema':
        [u'http://www.tardis.edu.au/schemas/ngs/project']}}  # type: dict
        
        # instrument_model fixture
        self._instrument_model__schema = {u'pk': None, u'model':
        u'tardis_portal.parametername', u'fields': {u'name':
        u'instrument_model', u'data_type': 2, u'immutable': True,
        u'is_searchable': True, u'choices': u'', u'comparison_type': 1,
        u'full_name': u'Instrument model', u'units': u'', u'order': 9999,
        u'schema': [u'http://www.tardis.edu.au/schemas/ngs/project']}}  # type: dict
        
        # read_cycles fixture
        self._read_cycles__schema = {u'pk': None, u'model':
        u'tardis_portal.parametername', u'fields': {u'name': u'read_cycles',
        u'data_type': 2, u'immutable': True, u'is_searchable': False,
        u'choices': u'', u'comparison_type': 1, u'full_name': u'Number of'
        u'cycles in each read [index reads in (brackets)]', u'units': u'',
        u'order': 9999, u'schema':
        [u'http://www.tardis.edu.au/schemas/ngs/project']}}  # type: dict
        
        # chemistry fixture
        self._chemistry__schema = {u'pk': None, u'model':
        u'tardis_portal.parametername', u'fields': {u'name': u'chemistry',
        u'data_type': 2, u'immutable': True, u'is_searchable': False,
        u'choices': u'', u'comparison_type': 1, u'full_name': u'Terminator'
        u'chemistry', u'units': u'', u'order': 9999, u'schema':
        [u'http://www.tardis.edu.au/schemas/ngs/project']}}  # type: dict
        
        # operator_name fixture
        self._operator_name__schema = {u'pk': None, u'model':
        u'tardis_portal.parametername', u'fields': {u'name': u'operator_name',
        u'data_type': 2, u'immutable': True, u'is_searchable': False,
        u'choices': u'', u'comparison_type': 1, u'full_name': u'Instrument'
        u'operator', u'units': u'', u'order': 9999, u'schema':
        [u'http://www.tardis.edu.au/schemas/ngs/project']}}  # type: dict
        
        # rta_version fixture
        self._rta_version__schema = {u'pk': None, u'model':
        u'tardis_portal.parametername', u'fields': {u'name': u'rta_version',
        u'data_type': 2, u'immutable': True, u'is_searchable': False,
        u'choices': u'', u'comparison_type': 1, u'full_name': u'Illumina RTA'
        u'version', u'units': u'', u'order': 9999, u'schema':
        [u'http://www.tardis.edu.au/schemas/ngs/project']}}  # type: dict
        

        self._subtype = "demultiplexed-samples"  # type: unicode
        self._model = "tardis_portal.schema"  # type: unicode
        self._name = "Sequencing Project (Demultiplexed Sample Set)"  # type: unicode
        self._pk = None  # type: NoneType
        self._type = 1  # type: int
        self._hidden = False  # type: bool
        self._namespace = "http://www.tardis.edu.au/schemas/ngs/project"  # type: unicode
        self._immutable = True  # type: bool


class NucleotideRawReadsBase(Dataset):
    """


    :type run_id: unicode
    :type project_experiment: unicode
    :type run_experiment: unicode
    :type fastqc_dataset: unicode
    :type run_number: float
    :type flowcell_id: unicode
    :type instrument_id: unicode
    :type instrument_model: unicode
    :type read_cycles: unicode
    :type chemistry: unicode
    :type operator_name: unicode
    :type rta_version: unicode
    """

    def __init__(self):
        super(Dataset, self).__init__()
        # Run Unique ID
        self.run_id = None  # type: unicode
        
        # Project Experiment link
        self.project_experiment = None  # type: unicode
        
        # Run Experiment link
        self.run_experiment = None  # type: unicode
        
        # Associated FastQC reports
        self.fastqc_dataset = None  # type: unicode
        
        # Run number
        self.run_number = None  # type: float
        
        # Flowcell ID
        self.flowcell_id = None  # type: unicode
        
        # Instrument ID
        self.instrument_id = None  # type: unicode
        
        # Instrument model
        self.instrument_model = None  # type: unicode
        
        # Number of cycles in each read [index reads in (brackets)]
        self.read_cycles = None  # type: unicode
        
        # Terminator chemistry
        self.chemistry = None  # type: unicode
        
        # Instrument operator
        self.operator_name = None  # type: unicode
        
        # Illumina RTA version
        self.rta_version = None  # type: unicode
        

        # Dictionaries to allow reconstitution of the schema for each parameter

        # run_id fixture
        self._run_id__schema = {u'pk': None, u'model':
        u'tardis_portal.parametername', u'fields': {u'name': u'run_id',
        u'data_type': 2, u'immutable': True, u'is_searchable': True,
        u'choices': u'', u'comparison_type': 1, u'full_name': u'Run Unique'
        u'ID', u'units': u'', u'order': 9999, u'schema':
        [u'http://www.tardis.edu.au/schemas/ngs/project/raw_reads']}}  # type: dict
        
        # project_experiment fixture
        self._project_experiment__schema = {u'pk': None, u'model':
        u'tardis_portal.parametername', u'fields': {u'name':
        u'project_experiment', u'data_type': 4, u'immutable': True,
        u'is_searchable': True, u'choices': u'', u'comparison_type': 1,
        u'full_name': u'Project Experiment link', u'units': u'', u'order':
        9999, u'schema':
        [u'http://www.tardis.edu.au/schemas/ngs/project/raw_reads']}}  # type: dict
        
        # run_experiment fixture
        self._run_experiment__schema = {u'pk': None, u'model':
        u'tardis_portal.parametername', u'fields': {u'name':
        u'run_experiment', u'data_type': 4, u'immutable': True,
        u'is_searchable': True, u'choices': u'', u'comparison_type': 1,
        u'full_name': u'Run Experiment link', u'units': u'', u'order': 9999,
        u'schema':
        [u'http://www.tardis.edu.au/schemas/ngs/project/raw_reads']}}  # type: dict
        
        # fastqc_dataset fixture
        self._fastqc_dataset__schema = {u'pk': None, u'model':
        u'tardis_portal.parametername', u'fields': {u'name':
        u'fastqc_dataset', u'data_type': 4, u'immutable': False,
        u'is_searchable': True, u'choices': u'', u'comparison_type': 1,
        u'full_name': u'Associated FastQC reports', u'units': u'', u'order':
        9999, u'schema':
        [u'http://www.tardis.edu.au/schemas/ngs/project/raw_reads']}}  # type: dict
        
        # run_number fixture
        self._run_number__schema = {u'pk': None, u'model':
        u'tardis_portal.parametername', u'fields': {u'name': u'run_number',
        u'data_type': 1, u'immutable': True, u'is_searchable': True,
        u'choices': u'', u'comparison_type': 1, u'full_name': u'Run number',
        u'units': u'', u'order': 9999, u'schema':
        [u'http://www.tardis.edu.au/schemas/ngs/project/raw_reads']}}  # type: dict
        
        # flowcell_id fixture
        self._flowcell_id__schema = {u'pk': None, u'model':
        u'tardis_portal.parametername', u'fields': {u'name': u'flowcell_id',
        u'data_type': 2, u'immutable': True, u'is_searchable': True,
        u'choices': u'', u'comparison_type': 1, u'full_name': u'Flowcell ID',
        u'units': u'', u'order': 9999, u'schema':
        [u'http://www.tardis.edu.au/schemas/ngs/project/raw_reads']}}  # type: dict
        
        # instrument_id fixture
        self._instrument_id__schema = {u'pk': None, u'model':
        u'tardis_portal.parametername', u'fields': {u'name': u'instrument_id',
        u'data_type': 2, u'immutable': True, u'is_searchable': True,
        u'choices': u'', u'comparison_type': 1, u'full_name': u'Instrument'
        u'ID', u'units': u'', u'order': 9999, u'schema':
        [u'http://www.tardis.edu.au/schemas/ngs/project/raw_reads']}}  # type: dict
        
        # instrument_model fixture
        self._instrument_model__schema = {u'pk': None, u'model':
        u'tardis_portal.parametername', u'fields': {u'name':
        u'instrument_model', u'data_type': 2, u'immutable': True,
        u'is_searchable': True, u'choices': u'', u'comparison_type': 1,
        u'full_name': u'Instrument model', u'units': u'', u'order': 9999,
        u'schema':
        [u'http://www.tardis.edu.au/schemas/ngs/project/raw_reads']}}  # type: dict
        
        # read_cycles fixture
        self._read_cycles__schema = {u'pk': None, u'model':
        u'tardis_portal.parametername', u'fields': {u'name': u'read_cycles',
        u'data_type': 2, u'immutable': True, u'is_searchable': False,
        u'choices': u'', u'comparison_type': 1, u'full_name': u'Number of'
        u'cycles in each read [index reads in (brackets)]', u'units': u'',
        u'order': 9999, u'schema':
        [u'http://www.tardis.edu.au/schemas/ngs/project/raw_reads']}}  # type: dict
        
        # chemistry fixture
        self._chemistry__schema = {u'pk': None, u'model':
        u'tardis_portal.parametername', u'fields': {u'name': u'chemistry',
        u'data_type': 2, u'immutable': True, u'is_searchable': False,
        u'choices': u'', u'comparison_type': 1, u'full_name': u'Terminator'
        u'chemistry', u'units': u'', u'order': 9999, u'schema':
        [u'http://www.tardis.edu.au/schemas/ngs/project/raw_reads']}}  # type: dict
        
        # operator_name fixture
        self._operator_name__schema = {u'pk': None, u'model':
        u'tardis_portal.parametername', u'fields': {u'name': u'operator_name',
        u'data_type': 2, u'immutable': True, u'is_searchable': False,
        u'choices': u'', u'comparison_type': 1, u'full_name': u'Instrument'
        u'operator', u'units': u'', u'order': 9999, u'schema':
        [u'http://www.tardis.edu.au/schemas/ngs/project/raw_reads']}}  # type: dict
        
        # rta_version fixture
        self._rta_version__schema = {u'pk': None, u'model':
        u'tardis_portal.parametername', u'fields': {u'name': u'rta_version',
        u'data_type': 2, u'immutable': True, u'is_searchable': False,
        u'choices': u'', u'comparison_type': 1, u'full_name': u'Illumina RTA'
        u'version', u'units': u'', u'order': 9999, u'schema':
        [u'http://www.tardis.edu.au/schemas/ngs/project/raw_reads']}}  # type: dict
        

        self._subtype = "nucleotide-raw-reads"  # type: unicode
        self._model = "tardis_portal.schema"  # type: unicode
        self._name = "Nucleotide Sequencing Project Raw Reads"  # type: unicode
        self._pk = None  # type: NoneType
        self._type = 2  # type: int
        self._hidden = False  # type: bool
        self._namespace = "http://www.tardis.edu.au/schemas/ngs/project/raw_reads"  # type: unicode
        self._immutable = True  # type: bool


class FastqcReportsBase(Dataset):
    """


    :type run_id: unicode
    :type project: unicode
    :type raw_reads_dataset: unicode
    :type fastqc_version: unicode
    """

    def __init__(self):
        super(Dataset, self).__init__()
        # Run Unique ID
        self.run_id = None  # type: unicode
        
        # Project name
        self.project = None  # type: unicode
        
        # Raw reads project link
        self.raw_reads_dataset = None  # type: unicode
        
        # FastQC software version
        self.fastqc_version = None  # type: unicode
        

        # Dictionaries to allow reconstitution of the schema for each parameter

        # run_id fixture
        self._run_id__schema = {u'pk': None, u'model':
        u'tardis_portal.parametername', u'fields': {u'name': u'run_id',
        u'data_type': 2, u'immutable': True, u'is_searchable': True,
        u'choices': u'', u'comparison_type': 1, u'full_name': u'Run Unique'
        u'ID', u'units': u'', u'order': 9999, u'schema':
        [u'http://www.tardis.edu.au/schemas/ngs/project/fastqc']}}  # type: dict
        
        # project fixture
        self._project__schema = {u'pk': None, u'model':
        u'tardis_portal.parametername', u'fields': {u'name': u'project',
        u'data_type': 2, u'immutable': True, u'is_searchable': False,
        u'choices': u'', u'comparison_type': 1, u'full_name': u'Project name',
        u'units': u'', u'order': 9999, u'schema':
        [u'http://www.tardis.edu.au/schemas/ngs/project/fastqc']}}  # type: dict
        
        # raw_reads_dataset fixture
        self._raw_reads_dataset__schema = {u'pk': None, u'model':
        u'tardis_portal.parametername', u'fields': {u'name':
        u'raw_reads_dataset', u'data_type': 4, u'immutable': False,
        u'is_searchable': True, u'choices': u'', u'comparison_type': 1,
        u'full_name': u'Raw reads project link', u'units': u'', u'order':
        9999, u'schema':
        [u'http://www.tardis.edu.au/schemas/ngs/project/fastqc']}}  # type: dict
        
        # fastqc_version fixture
        self._fastqc_version__schema = {u'pk': None, u'model':
        u'tardis_portal.parametername', u'fields': {u'name':
        u'fastqc_version', u'data_type': 2, u'immutable': True,
        u'is_searchable': False, u'choices': u'', u'comparison_type': 1,
        u'full_name': u'FastQC software version', u'units': u'', u'order':
        9999, u'schema':
        [u'http://www.tardis.edu.au/schemas/ngs/project/fastqc']}}  # type: dict
        

        self._subtype = "fastqc-reports"  # type: unicode
        self._model = "tardis_portal.schema"  # type: unicode
        self._name = "FastQC Reports"  # type: unicode
        self._pk = None  # type: NoneType
        self._type = 2  # type: int
        self._hidden = False  # type: bool
        self._namespace = "http://www.tardis.edu.au/schemas/ngs/project/fastqc"  # type: unicode
        self._immutable = True  # type: bool


class HiddenFastqcProjectSummaryBase(Dataset):
    """


    :type _hidden_fastqc_summary_json: dict
    :type fastqc_version: unicode
    """

    def __init__(self):
        super(Dataset, self).__init__()
        # (Hidden) FastQC summary for all samples (JSON)
        self._hidden_fastqc_summary_json = None  # type: dict
        
        # FastQC software version
        self.fastqc_version = None  # type: unicode
        

        # Dictionaries to allow reconstitution of the schema for each parameter

        # _hidden_fastqc_summary_json fixture
        self.__hidden_fastqc_summary_json__schema = {u'pk': None, u'model':
        u'tardis_portal.parametername', u'fields': {u'name':
        u'__hidden__fastqc_summary_json', u'data_type': 8, u'immutable': True,
        u'is_searchable': False, u'choices': u'', u'comparison_type': 1,
        u'full_name': u'(Hidden) FastQC summary for all samples (JSON)',
        u'units': u'fastqc-summary-table', u'order': 9999, u'schema': [u'http:'
        u'//www.tardis.edu.au/schemas/ngs/project/__hidden__fastqc_summary']}}  # type: dict
        
        # fastqc_version fixture
        self._fastqc_version__schema = {u'pk': None, u'model':
        u'tardis_portal.parametername', u'fields': {u'name':
        u'fastqc_version', u'data_type': 2, u'immutable': True,
        u'is_searchable': False, u'choices': u'', u'comparison_type': 1,
        u'full_name': u'FastQC software version', u'units': u'', u'order':
        9999, u'schema': [u'http://www.tardis.edu.au/schemas/ngs/project/__hid'
        u'den__fastqc_summary']}}  # type: dict
        

        self._subtype = "hidden-fastqc-project-summary"  # type: unicode
        self._model = "tardis_portal.schema"  # type: unicode
        self._name = "FastQC Project Summary"  # type: unicode
        self._pk = None  # type: NoneType
        self._type = 2  # type: int
        self._hidden = True  # type: bool
        self._namespace = "http://www.tardis.edu.au/schemas/ngs/project/__hidden__fastqc_summary"  # type: unicode
        self._immutable = True  # type: bool


class FastqRawReadsBase(DataFile):
    """


    :type run_id: unicode
    :type sample_id: unicode
    :type reference_genome: unicode
    :type index_sequence: unicode
    :type is_control: unicode
    :type recipe: unicode
    :type operator_name: unicode
    :type description: unicode
    :type project: unicode
    :type number_of_reads: float
    :type number_of_poor_quality_reads: float
    :type read_length: float
    """

    def __init__(self):
        super(DataFile, self).__init__()
        # Run Unique ID
        self.run_id = None  # type: unicode
        
        # Sample ID
        self.sample_id = None  # type: unicode
        
        # Reference genome
        self.reference_genome = None  # type: unicode
        
        # Index sequence (barcode) for this sample
        self.index_sequence = None  # type: unicode
        
        # Is control ?
        self.is_control = None  # type: unicode
        
        # Recipe
        self.recipe = None  # type: unicode
        
        # Instrument Operator
        self.operator_name = None  # type: unicode
        
        # Description
        self.description = None  # type: unicode
        
        # Project name
        self.project = None  # type: unicode
        
        # Number of reads
        self.number_of_reads = None  # type: float
        
        # Number of reads flagged as poor quality (FastQC)
        self.number_of_poor_quality_reads = None  # type: float
        
        # Read length
        self.read_length = None  # type: float
        

        # Dictionaries to allow reconstitution of the schema for each parameter

        # run_id fixture
        self._run_id__schema = {u'pk': None, u'model':
        u'tardis_portal.parametername', u'fields': {u'name': u'run_id',
        u'data_type': 2, u'immutable': True, u'is_searchable': True,
        u'choices': u'', u'comparison_type': 1, u'full_name': u'Run Unique'
        u'ID', u'units': u'', u'order': 9999, u'schema':
        [u'http://www.tardis.edu.au/schemas/ngs/file/fastq']}}  # type: dict
        
        # sample_id fixture
        self._sample_id__schema = {u'pk': None, u'model':
        u'tardis_portal.parametername', u'fields': {u'name': u'sample_id',
        u'data_type': 2, u'immutable': True, u'is_searchable': False,
        u'choices': u'', u'comparison_type': 1, u'full_name': u'Sample ID',
        u'units': u'', u'order': 9999, u'schema':
        [u'http://www.tardis.edu.au/schemas/ngs/file/fastq']}}  # type: dict
        
        # reference_genome fixture
        self._reference_genome__schema = {u'pk': None, u'model':
        u'tardis_portal.parametername', u'fields': {u'name':
        u'reference_genome', u'data_type': 2, u'immutable': True,
        u'is_searchable': False, u'choices': u'', u'comparison_type': 1,
        u'full_name': u'Reference genome', u'units': u'', u'order': 9999,
        u'schema': [u'http://www.tardis.edu.au/schemas/ngs/file/fastq']}}  # type: dict
        
        # index_sequence fixture
        self._index_sequence__schema = {u'pk': None, u'model':
        u'tardis_portal.parametername', u'fields': {u'name':
        u'index_sequence', u'data_type': 2, u'immutable': True,
        u'is_searchable': False, u'choices': u'', u'comparison_type': 1,
        u'full_name': u'Index sequence (barcode) for this sample', u'units':
        u'', u'order': 9999, u'schema':
        [u'http://www.tardis.edu.au/schemas/ngs/file/fastq']}}  # type: dict
        
        # is_control fixture
        self._is_control__schema = {u'pk': None, u'model':
        u'tardis_portal.parametername', u'fields': {u'name': u'is_control',
        u'data_type': 2, u'immutable': True, u'is_searchable': False,
        u'choices': u'', u'comparison_type': 1, u'full_name': u'Is control ?',
        u'units': u'', u'order': 9999, u'schema':
        [u'http://www.tardis.edu.au/schemas/ngs/file/fastq']}}  # type: dict
        
        # recipe fixture
        self._recipe__schema = {u'pk': None, u'model':
        u'tardis_portal.parametername', u'fields': {u'name': u'recipe',
        u'data_type': 2, u'immutable': True, u'is_searchable': False,
        u'choices': u'', u'comparison_type': 1, u'full_name': u'Recipe',
        u'units': u'', u'order': 9999, u'schema':
        [u'http://www.tardis.edu.au/schemas/ngs/file/fastq']}}  # type: dict
        
        # operator_name fixture
        self._operator_name__schema = {u'pk': None, u'model':
        u'tardis_portal.parametername', u'fields': {u'name': u'operator_name',
        u'data_type': 2, u'immutable': True, u'is_searchable': False,
        u'choices': u'', u'comparison_type': 1, u'full_name': u'Instrument'
        u'Operator', u'units': u'', u'order': 9999, u'schema':
        [u'http://www.tardis.edu.au/schemas/ngs/file/fastq']}}  # type: dict
        
        # description fixture
        self._description__schema = {u'pk': None, u'model':
        u'tardis_portal.parametername', u'fields': {u'name': u'description',
        u'data_type': 2, u'immutable': True, u'is_searchable': False,
        u'choices': u'', u'comparison_type': 1, u'full_name': u'Description',
        u'units': u'', u'order': 9999, u'schema':
        [u'http://www.tardis.edu.au/schemas/ngs/file/fastq']}}  # type: dict
        
        # project fixture
        self._project__schema = {u'pk': None, u'model':
        u'tardis_portal.parametername', u'fields': {u'name': u'project',
        u'data_type': 2, u'immutable': True, u'is_searchable': False,
        u'choices': u'', u'comparison_type': 1, u'full_name': u'Project name',
        u'units': u'', u'order': 9999, u'schema':
        [u'http://www.tardis.edu.au/schemas/ngs/file/fastq']}}  # type: dict
        
        # number_of_reads fixture
        self._number_of_reads__schema = {u'pk': None, u'model':
        u'tardis_portal.parametername', u'fields': {u'name':
        u'number_of_reads', u'data_type': 1, u'immutable': True,
        u'is_searchable': False, u'choices': u'', u'comparison_type': 1,
        u'full_name': u'Number of reads', u'units': u'', u'order': 9999,
        u'schema': [u'http://www.tardis.edu.au/schemas/ngs/file/fastq']}}  # type: dict
        
        # number_of_poor_quality_reads fixture
        self._number_of_poor_quality_reads__schema = {u'pk': None, u'model':
        u'tardis_portal.parametername', u'fields': {u'name':
        u'number_of_poor_quality_reads', u'data_type': 1, u'immutable': True,
        u'is_searchable': False, u'choices': u'', u'comparison_type': 1,
        u'full_name': u'Number of reads flagged as poor quality (FastQC)',
        u'units': u'', u'order': 9999, u'schema':
        [u'http://www.tardis.edu.au/schemas/ngs/file/fastq']}}  # type: dict
        
        # read_length fixture
        self._read_length__schema = {u'pk': None, u'model':
        u'tardis_portal.parametername', u'fields': {u'name': u'read_length',
        u'data_type': 1, u'immutable': True, u'is_searchable': False,
        u'choices': u'', u'comparison_type': 1, u'full_name': u'Read length',
        u'units': u'', u'order': 9999, u'schema':
        [u'http://www.tardis.edu.au/schemas/ngs/file/fastq']}}  # type: dict
        

        self._subtype = "fastq-raw-reads"  # type: unicode
        self._model = "tardis_portal.schema"  # type: unicode
        self._name = "Nucleotide Sequence Raw Reads (FASTQ)"  # type: unicode
        self._pk = None  # type: NoneType
        self._type = 3  # type: int
        self._hidden = False  # type: bool
        self._namespace = "http://www.tardis.edu.au/schemas/ngs/file/fastq"  # type: unicode
        self._immutable = True  # type: bool


class FastqcOutputBase(DataFile):
    """


    :type run_id: unicode
    :type sample_id: unicode
    :type project: unicode
    :type fastqc_version: unicode
    """

    def __init__(self):
        super(DataFile, self).__init__()
        # Run Unique ID
        self.run_id = None  # type: unicode
        
        # Sample ID
        self.sample_id = None  # type: unicode
        
        # Project name
        self.project = None  # type: unicode
        
        # FastQC software version
        self.fastqc_version = None  # type: unicode
        

        # Dictionaries to allow reconstitution of the schema for each parameter

        # run_id fixture
        self._run_id__schema = {u'pk': None, u'model':
        u'tardis_portal.parametername', u'fields': {u'name': u'run_id',
        u'data_type': 2, u'immutable': True, u'is_searchable': True,
        u'choices': u'', u'comparison_type': 1, u'full_name': u'Run Unique'
        u'ID', u'units': u'', u'order': 9999, u'schema':
        [u'http://www.tardis.edu.au/schemas/ngs/file/fastqc']}}  # type: dict
        
        # sample_id fixture
        self._sample_id__schema = {u'pk': None, u'model':
        u'tardis_portal.parametername', u'fields': {u'name': u'sample_id',
        u'data_type': 2, u'immutable': True, u'is_searchable': False,
        u'choices': u'', u'comparison_type': 1, u'full_name': u'Sample ID',
        u'units': u'', u'order': 9999, u'schema':
        [u'http://www.tardis.edu.au/schemas/ngs/file/fastqc']}}  # type: dict
        
        # project fixture
        self._project__schema = {u'pk': None, u'model':
        u'tardis_portal.parametername', u'fields': {u'name': u'project',
        u'data_type': 2, u'immutable': True, u'is_searchable': False,
        u'choices': u'', u'comparison_type': 1, u'full_name': u'Project name',
        u'units': u'', u'order': 9999, u'schema':
        [u'http://www.tardis.edu.au/schemas/ngs/file/fastqc']}}  # type: dict
        
        # fastqc_version fixture
        self._fastqc_version__schema = {u'pk': None, u'model':
        u'tardis_portal.parametername', u'fields': {u'name':
        u'fastqc_version', u'data_type': 2, u'immutable': True,
        u'is_searchable': False, u'choices': u'', u'comparison_type': 1,
        u'full_name': u'FastQC software version', u'units': u'', u'order':
        9999, u'schema':
        [u'http://www.tardis.edu.au/schemas/ngs/file/fastqc']}}  # type: dict
        

        self._subtype = "fastqc-output"  # type: unicode
        self._model = "tardis_portal.schema"  # type: unicode
        self._name = "FastQC report"  # type: unicode
        self._pk = None  # type: NoneType
        self._type = 3  # type: int
        self._hidden = False  # type: bool
        self._namespace = "http://www.tardis.edu.au/schemas/ngs/file/fastqc"  # type: unicode
        self._immutable = True  # type: bool

