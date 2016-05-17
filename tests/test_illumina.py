import os
import unittest
from datetime import datetime
from semantic_version import Version as SemanticVersion
import mytardis_ngs_ingestor
from mytardis_ngs_ingestor.illumina.run_info import \
    parse_samplesheet, \
    filter_samplesheet_by_project, \
    filter_samplesheet_by_project, \
    rta_complete_parser


class IlluminaParserTestCase(unittest.TestCase):
    def setUp(self):
        self.run1_dir = os.path.join(
            os.path.dirname(__file__),
            'test_data/runs/130907_SNL177_0001_AH9PJLADXZ')

        self.run2_dir = os.path.join(
            os.path.dirname(__file__),
            'test_data/runs/140907_SNL177_0002_AC9246ACXX')

        self.samplesheet_csv_path = os.path.join(os.path.dirname(__file__),
                                                 self.run1_dir,
                                                 'SampleSheet.csv')
        self.samplesheet_v4_path = os.path.join(os.path.dirname(__file__),
                                                self.run2_dir,
                                                'SampleSheet.csv')

    def tearDown(self):
        pass

    def test_parse_samplesheet_v4(self):
        samples, chemistry = parse_samplesheet(self.samplesheet_v4_path,
                                               standardize_keys=True)
        self.assertEqual(chemistry, 'TruSeq LT')

        expected = [
            {'index': 'CGATGT', 'Lane': '1', 'Description': '',
             'SampleID': '16-00982', 'SamplePlate': '', 'I7IndexID': 'A002',
             'SampleWell': '', 'SampleProject': 'StephenLavelle',
             'SampleName': 'QQ1H2O1'},
            {'index': 'TGACCA', 'Lane': '1', 'Description': '',
             'SampleID': '16-00983', 'SamplePlate': '', 'I7IndexID': 'A004',
             'SampleWell': '', 'SampleProject': 'StephenLavelle',
             'SampleName': 'QQ1H2O2'},
            {'index': 'CAGATC', 'Lane': '1', 'Description': '',
             'SampleID': '16-00984', 'SamplePlate': '', 'I7IndexID': 'A007',
             'SampleWell': '', 'SampleProject': 'StephenLavelle',
             'SampleName': 'QQ1H2O3'},
            {'index': 'AACCAG', 'Lane': '2', 'Description': '',
             'SampleID': '16-00487', 'SamplePlate': '', 'I7IndexID': 'A001',
             'SampleWell': '', 'SampleProject': 'Shigeru_Miyamoto',
             'SampleName': 'QQInputF2'},
            {'index': 'TGGTGA', 'Lane': '2', 'Description': '',
             'SampleID': '16-00488', 'SamplePlate': '', 'I7IndexID': 'A002',
             'SampleWell': '', 'SampleProject': 'Shigeru_Miyamoto',
             'SampleName': 'QQH4K4F2'},
            {'index': 'AGTGAG', 'Lane': '2', 'Description': '',
             'SampleID': '16-00489', 'SamplePlate': '', 'I7IndexID': 'A003',
             'SampleWell': '', 'SampleProject': 'Shigeru_Miyamoto',
             'SampleName': 'QQH4K9F2'},
            {'index': 'AACCAG', 'Lane': '3', 'Description': '',
             'SampleID': '16-01787', 'SamplePlate': '', 'I7IndexID': 'A001',
             'SampleWell': '', 'SampleProject': 'Phr00t',
             'SampleName': 'Q1N'},
            {'index': 'TGGTGA', 'Lane': '3', 'Description': '',
             'SampleID': '16-01788', 'SamplePlate': '', 'I7IndexID': 'A002',
             'SampleWell': '', 'SampleProject': 'Phr00t',
             'SampleName': 'Q1L'},
            {'index': 'AACCAG', 'Lane': '4', 'Description': '',
             'SampleID': '16-01787', 'SamplePlate': '', 'I7IndexID': 'A001',
             'SampleWell': '', 'SampleProject': 'Phr00t',
             'SampleName': 'Q1N'}
        ]

        for expected_sample, sample in zip(expected, samples):
            self.assertDictEqual(sample, expected_sample)

        self.assertEqual(chemistry, 'TruSeq LT')

    def test_parse_samplesheet_csv(self):
        samples, chemistry = parse_samplesheet(self.samplesheet_csv_path)
        expected = [
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
             'SampleProject': 'Walter_White', 'Recipe': '',
             'SampleID': '14-05655-SW38', 'FCID': 'H9PJLADXZ',
             'SampleRef': 'Hg19', 'Operator': 'TW'},
            {'Control': 'N', 'Index': 'AAGAGG', 'Lane': '2', 'Description': '',
             'SampleProject': 'Walter_White', 'Recipe': '',
             'SampleID': '14-05658-SW41', 'FCID': 'H9PJLADXZ',
             'SampleRef': 'Hg19', 'Operator': 'TW'},
            {'Control': 'N', 'Index': 'GGAGAA', 'Lane': '2', 'Description': '',
             'SampleProject': 'Walter_White', 'Recipe': '',
             'SampleID': '14-05659-SW42', 'FCID': 'H9PJLADXZ',
             'SampleRef': 'Hg19', 'Operator': 'TW'},
            {'Control': 'N', 'Index': 'AGCATG', 'Lane': '2', 'Description': '',
             'SampleProject': 'Walter_White', 'Recipe': '',
             'SampleID': '14-05660-SW43', 'FCID': 'H9PJLADXZ',
             'SampleRef': 'Hg19', 'Operator': 'TW'},
            {'Control': 'N', 'Index': 'GAGTCA', 'Lane': '2', 'Description': '',
             'SampleProject': 'Walter_White', 'Recipe': '',
             'SampleID': '14-05661-SW44', 'FCID': 'H9PJLADXZ',
             'SampleRef': 'Hg19', 'Operator': 'TW'},
            {'Control': 'N', 'Index': 'CGTAGA', 'Lane': '2', 'Description': '',
             'SampleProject': 'Walter_White', 'Recipe': '',
             'SampleID': '14-06203-SW45', 'FCID': 'H9PJLADXZ',
             'SampleRef': 'Hg19', 'Operator': 'TW'},
            {'Control': 'N', 'Index': 'TCAGAG', 'Lane': '2', 'Description': '',
             'SampleProject': 'Walter_White', 'Recipe': '',
             'SampleID': '14-06204-SW46', 'FCID': 'H9PJLADXZ',
             'SampleRef': 'Hg19', 'Operator': 'TW'}]

        for sample_line, expected_line in zip(samples, expected):
            self.assertDictEqual(sample_line, expected_line)

        self.assertEqual(chemistry, 'TruSeq LT')

    def test_filter_samplesheet_by_project(self):
        project_lines = filter_samplesheet_by_project(
            self.samplesheet_csv_path, 'GusFring')
        expected = [
            'FCID,Lane,SampleID,SampleRef,Index,Description,Control,Recipe,Operator,SampleProject\r\n',
            'H9PJLADXZ,1,14-06205-OCT4-5,Hg19,AACCAG,May contain nuts,N,,TW,GusFring\r\n',
            'H9PJLADXZ,1,14-06206-OCT4-15,Hg19,TGGTGA,,N,,TW,GusFring\r\n',
            'H9PJLADXZ,1,14-06207-ZAX-5,Hg19,AGTGAG,,N,,TW,GusFring\r\n',
            'H9PJLADXZ,1,14-06208-ZAX-15,Hg19,GCACTA,,N,,TW,GusFring\r\n',
            'H9PJLADXZ,1,14-06200-Input,Hg19,TTGGCA,,N,,TW,GusFring\r\n',
            '#_IEMVERSION_3_TruSeq LT,,,,,,,,,\r\n'
        ]

        for expected_line, project_line in zip(expected, project_lines):
            self.assertEqual(project_line, expected_line)

    # TODO: These times and RTA versions aren't atually consistent
    #       with the other times of the mock runs. Make them
    #       consistent.
    def test_rta_complete_parser(self):
        # 1.x = '9/7/2013,18:12:53.149,Illumina RTA 1.18.64'
        date, version = rta_complete_parser(self.run1_dir)
        self.assertEqual(date, datetime(2013, 9, 7, 18, 12, 53, 149000))
        self.assertEqual(version, 'Illumina RTA 1.17.20.0')

        # 2.x = 'RTA 2.7.3 completed on 9/7/2014 3:31:22 AM'
        date, version = rta_complete_parser(self.run2_dir)
        self.assertEqual(date, datetime(2014, 9, 7, 3, 31, 22))
        self.assertEqual(version, 'RTA 2.7.3')


class VersionTest(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_semver(self):
        # parse the version number to ensure it's valid (invalid version number
        # formats will raise a ValueError)
        ingestor_version = SemanticVersion(mytardis_ngs_ingestor.__version__)
        self.assertEqual(str(ingestor_version),
                         mytardis_ngs_ingestor.__version__)


if __name__ == '__main__':
    unittest.main()
