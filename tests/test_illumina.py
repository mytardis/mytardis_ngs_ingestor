import os
import unittest
from semantic_version import Version as SemanticVersion
import mytardis_uploader
from mytardis_uploader.illumina_uploader import \
    parse_samplesheet, \
    filter_samplesheet_by_project, \
    filter_samplesheet_by_project


class IlluminaParserTestCase(unittest.TestCase):
    def setUp(self):
        self.samplesheet_path = os.path.join(os.path.dirname(__file__),
                                             os.path.join('test_data',
                                                          'SampleSheet.csv'))

    def tearDown(self):
        pass

    def test_parse_samplesheet(self):
        samples, chemistry = parse_samplesheet(self.samplesheet_path)
        self.assertEqual(samples, [
            {'Control': 'N', 'Index': 'AACCAG', 'Lane': '1',
             'Description': 'May contain nuts',
             'SampleProject': 'GusFring', 'Recipe': '',
             'SampleID': '14-06205-OCT4-5', 'FCID': 'H9PJLADXZ',
             'SampleRef': 'Hg19', 'Operator': 'TW'},
            {'Control': 'N', 'Index': 'TGGTGA', 'Lane': '1', 'Description': '',
             'SampleProject': 'GusFring', 'Recipe': '',
             'SampleID': '14-06206-OCT4-15', 'FCID': 'H9PJLADXZ',
             'SampleRef': 'Hg19', 'Operator': 'TW'},
            {'Control': 'N', 'Index': 'AGTGAG', 'Lane': '1', 'Description': '',
             'SampleProject': 'GusFring', 'Recipe': '',
             'SampleID': '14-06207-ZAX-5', 'FCID': 'H9PJLADXZ',
             'SampleRef': 'Hg19', 'Operator': 'TW'},
            {'Control': 'N', 'Index': 'GCACTA', 'Lane': '1', 'Description': '',
             'SampleProject': 'GusFring', 'Recipe': '',
             'SampleID': '14-06208-ZAX-15', 'FCID': 'H9PJLADXZ',
             'SampleRef': 'Hg19', 'Operator': 'TW'},
            {'Control': 'N', 'Index': 'TTGGCA', 'Lane': '1', 'Description': '',
             'SampleProject': 'GusFring', 'Recipe': '',
             'SampleID': '14-06200-Input', 'FCID': 'H9PJLADXZ',
             'SampleRef': 'Hg19', 'Operator': 'TW'},
            {'Control': 'N', 'Index': 'ACCTCA', 'Lane': '2', 'Description': '',
             'SampleProject': 'WalterWhite', 'Recipe': '',
             'SampleID': '14-05655-SW38', 'FCID': 'H9PJLADXZ',
             'SampleRef': 'Hg19', 'Operator': 'TW'},
            {'Control': 'N', 'Index': 'AAGAGG', 'Lane': '2', 'Description': '',
             'SampleProject': 'WalterWhite', 'Recipe': '',
             'SampleID': '14-05658-SW41', 'FCID': 'H9PJLADXZ',
             'SampleRef': 'Hg19', 'Operator': 'TW'},
            {'Control': 'N', 'Index': 'GGAGAA', 'Lane': '2', 'Description': '',
             'SampleProject': 'WalterWhite', 'Recipe': '',
             'SampleID': '14-05659-SW42', 'FCID': 'H9PJLADXZ',
             'SampleRef': 'Hg19', 'Operator': 'TW'},
            {'Control': 'N', 'Index': 'AGCATG', 'Lane': '2', 'Description': '',
             'SampleProject': 'WalterWhite', 'Recipe': '',
             'SampleID': '14-05660-SW43', 'FCID': 'H9PJLADXZ',
             'SampleRef': 'Hg19', 'Operator': 'TW'},
            {'Control': 'N', 'Index': 'GAGTCA', 'Lane': '2', 'Description': '',
             'SampleProject': 'WalterWhite', 'Recipe': '',
             'SampleID': '14-05661-SW44', 'FCID': 'H9PJLADXZ',
             'SampleRef': 'Hg19', 'Operator': 'TW'},
            {'Control': 'N', 'Index': 'CGTAGA', 'Lane': '2', 'Description': '',
             'SampleProject': 'WalterWhite', 'Recipe': '',
             'SampleID': '14-06203-SW45', 'FCID': 'H9PJLADXZ',
             'SampleRef': 'Hg19', 'Operator': 'TW'},
            {'Control': 'N', 'Index': 'TCAGAG', 'Lane': '2', 'Description': '',
             'SampleProject': 'WalterWhite', 'Recipe': '',
             'SampleID': '14-06204-SW46', 'FCID': 'H9PJLADXZ',
             'SampleRef': 'Hg19', 'Operator': 'TW'}])
        self.assertEqual(chemistry, 'TruSeq LT')

    def test_filter_samplesheet_by_project(self):
        project_lines = filter_samplesheet_by_project(
            self.samplesheet_path, 'GusFring')
        self.assertEqual(project_lines, [
            'FCID,Lane,SampleID,SampleRef,Index,Description,Control,Recipe,Operator,SampleProject\r\n',
            'H9PJLADXZ,1,14-06205-OCT4-5,Hg19,AACCAG,May contain nuts,N,,TW,GusFring\n',
            'H9PJLADXZ,1,14-06206-OCT4-15,Hg19,TGGTGA,,N,,TW,GusFring\n',
            'H9PJLADXZ,1,14-06207-ZAX-5,Hg19,AGTGAG,,N,,TW,GusFring\n',
            'H9PJLADXZ,1,14-06208-ZAX-15,Hg19,GCACTA,,N,,TW,GusFring\n',
            'H9PJLADXZ,1,14-06200-Input,Hg19,TTGGCA,,N,,TW,GusFring\n',
            '#_IEMVERSION_3_TruSeq LT,,,,,,,,,'])


class VersionTest(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_semver(self):
        # parse the version number to ensure it's valid (invalid version number
        # formats will raise a ValueError)
        ingestor_version = SemanticVersion(mytardis_uploader.__version__)
        self.assertEqual(str(ingestor_version), mytardis_uploader.__version__)


if __name__ == '__main__':
    unittest.main()
