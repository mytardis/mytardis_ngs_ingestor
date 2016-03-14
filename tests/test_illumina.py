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
        self.samplesheet_csv_path = os.path.join(os.path.dirname(__file__),
                                                 os.path.join(
                                                     'test_data',
                                                     'SampleSheet.csv'))
        self.samplesheet_v4_path = os.path.join(os.path.dirname(__file__),
                                                os.path.join(
                                                    'test_data',
                                                    'SampleSheet_v4.csv'))

    def tearDown(self):
        pass

    def test_parse_samplesheet_v4(self):
        samples, chemistry = parse_samplesheet(self.samplesheet_v4_path,
                                               standardize_keys=True)
        self.assertEqual(chemistry, 'TruSeq LT')
        self.assertEqual(samples, [
            {'SampleWell': '', 'index': 'CGATGT', 'Lane': '1',
             'SamplePlate': '', 'I7IndexID': 'A002', 'SampleName': 'QQ1H2O1',
             'SampleProject': 'StephenLavelle', 'SampleID': '16-00982',
             'Description': ''},
            {'SampleWell': '', 'index': 'TGACCA', 'Lane': '1',
             'SamplePlate': '', 'I7IndexID': 'A004', 'SampleName': 'QQ1H2O2',
             'SampleProject': 'StephenLavelle', 'SampleID': '16-00983',
             'Description': ''},
            {'SampleWell': '', 'index': 'CAGATC', 'Lane': '1',
             'SamplePlate': '', 'I7IndexID': 'A007', 'SampleName': 'QQ1H2O3',
             'SampleProject': 'StephenLavelle', 'SampleID': '16-00984',
             'Description': ''},
            {'SampleWell': '', 'index': 'AACCAG', 'Lane': '5',
             'SamplePlate': '', 'I7IndexID': 'A001', 'SampleName': 'QQInputF2',
             'SampleProject': 'ShigeruMiyamoto', 'SampleID': '16-00487',
             'Description': ''},
            {'SampleWell': '', 'index': 'TGGTGA', 'Lane': '5',
             'SamplePlate': '', 'I7IndexID': 'A002', 'SampleName': 'QQH4K4F2',
             'SampleProject': 'ShigeruMiyamoto', 'SampleID': '16-00488',
             'Description': ''},
            {'SampleWell': '', 'index': 'AGTGAG', 'Lane': '5',
             'SamplePlate': '', 'I7IndexID': 'A003', 'SampleName': 'QQH4K9F2',
             'SampleProject': 'ShigeruMiyamoto', 'SampleID': '16-00489',
             'Description': ''},
            {'SampleWell': '', 'index': 'AACCAG', 'Lane': '7',
             'SamplePlate': '', 'I7IndexID': 'A001', 'SampleName': 'Q1N',
             'SampleProject': 'Phr00t', 'SampleID': '16-01787',
             'Description': ''},
            {'SampleWell': '', 'index': 'TGGTGA', 'Lane': '7',
             'SamplePlate': '', 'I7IndexID': 'A002', 'SampleName': 'Q1L',
             'SampleProject': 'Phr00t', 'SampleID': '16-01788',
             'Description': ''},
            {'SampleWell': '', 'index': 'AGTGAG', 'Lane': '7',
             'SamplePlate': '', 'I7IndexID': 'A003', 'SampleName': 'Q1H',
             'SampleProject': 'Phr00t', 'SampleID': '16-01789',
             'Description': ''}])

        samples, chemistry = parse_samplesheet(self.samplesheet_v4_path,
                                               standardize_keys=False)

        self.assertEqual(chemistry, 'TruSeq LT')
        self.assertEqual(samples, [
            {'index': 'CGATGT', 'Lane': '1', 'Description': '',
             'Sample_ID': '16-00982', 'Sample_Plate': '', 'I7_Index_ID': 'A002',
             'Sample_Well': '', 'Sample_Project': 'StephenLavelle',
             'Sample_Name': 'QQ1H2O1'},
            {'index': 'TGACCA', 'Lane': '1', 'Description': '',
             'Sample_ID': '16-00983', 'Sample_Plate': '', 'I7_Index_ID': 'A004',
             'Sample_Well': '', 'Sample_Project': 'StephenLavelle',
             'Sample_Name': 'QQ1H2O2'},
            {'index': 'CAGATC', 'Lane': '1', 'Description': '',
             'Sample_ID': '16-00984', 'Sample_Plate': '', 'I7_Index_ID': 'A007',
             'Sample_Well': '', 'Sample_Project': 'StephenLavelle',
             'Sample_Name': 'QQ1H2O3'},
            {'index': 'AACCAG', 'Lane': '5', 'Description': '',
             'Sample_ID': '16-00487', 'Sample_Plate': '', 'I7_Index_ID': 'A001',
             'Sample_Well': '', 'Sample_Project': 'ShigeruMiyamoto',
             'Sample_Name': 'QQInputF2'},
            {'index': 'TGGTGA', 'Lane': '5', 'Description': '',
             'Sample_ID': '16-00488', 'Sample_Plate': '', 'I7_Index_ID': 'A002',
             'Sample_Well': '', 'Sample_Project': 'ShigeruMiyamoto',
             'Sample_Name': 'QQH4K4F2'},
            {'index': 'AGTGAG', 'Lane': '5', 'Description': '',
             'Sample_ID': '16-00489', 'Sample_Plate': '', 'I7_Index_ID': 'A003',
             'Sample_Well': '', 'Sample_Project': 'ShigeruMiyamoto',
             'Sample_Name': 'QQH4K9F2'},
            {'index': 'AACCAG', 'Lane': '7', 'Description': '',
             'Sample_ID': '16-01787', 'Sample_Plate': '', 'I7_Index_ID': 'A001',
             'Sample_Well': '', 'Sample_Project': 'Phr00t',
             'Sample_Name': 'Q1N'},
            {'index': 'TGGTGA', 'Lane': '7', 'Description': '',
             'Sample_ID': '16-01788', 'Sample_Plate': '', 'I7_Index_ID': 'A002',
             'Sample_Well': '', 'Sample_Project': 'Phr00t',
             'Sample_Name': 'Q1L'},
            {'index': 'AGTGAG', 'Lane': '7', 'Description': '',
             'Sample_ID': '16-01789', 'Sample_Plate': '', 'I7_Index_ID': 'A003',
             'Sample_Well': '', 'Sample_Project': 'Phr00t',
             'Sample_Name': 'Q1H'}])

    def test_parse_samplesheet_csv(self):
        samples, chemistry = parse_samplesheet(self.samplesheet_csv_path)
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
            self.samplesheet_csv_path, 'GusFring')
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
