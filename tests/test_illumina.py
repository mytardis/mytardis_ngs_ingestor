import os
from os import path
import unittest
from datetime import datetime
from collections import OrderedDict
from semantic_version import Version as SemanticVersion
import mytardis_ngs_ingestor
from mytardis_ngs_ingestor.illumina.run_info import \
    parse_samplesheet, \
    filter_samplesheet_by_project, \
    filter_samplesheet_by_project, \
    rta_complete_parser, get_sample_project_mapping


class IlluminaParserTestCase(unittest.TestCase):
    def setUp(self):
        self.run1_dir = path.join(
            path.dirname(__file__),
            'test_data/runs/130907_DMO177_0001_AH9PJLADXZ')

        self.run2_dir = path.join(
            path.dirname(__file__),
            'test_data/runs/140907_DMO177_0002_AC9246ACXX')

        self.run3_dir = path.join(
            path.dirname(__file__),
            'test_data/runs/150907_M04242_0003_000000000-ANV1L')

        self.samplesheet_csv_path = path.join(os.path.dirname(__file__),
                                              self.run1_dir,
                                              'SampleSheet.csv')
        self.samplesheet_v4_path = path.join(os.path.dirname(__file__),
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

    def test_get_sample_project_mapping(self):
        bcl2fastq_output_path = path.join(self.run2_dir,
                                          'Data/Intensities/BaseCalls')
        mapping = get_sample_project_mapping(bcl2fastq_output_path)
        expected = OrderedDict(
            [(u'Phr00t', [u'Phr00t/16-01787/Q1N_S7_L003_R1_001.fastq.gz',
                          u'Phr00t/16-01787/Q1N_S7_L004_R1_001.fastq.gz',
                          u'Phr00t/16-01788/Q1L_S8_L003_R1_001.fastq.gz']),
             (u'Shigeru_Miyamoto', [
                 u'Shigeru_Miyamoto/16-00487/QQInputF2_S4_L002_R1_001.fastq.gz',
                 u'Shigeru_Miyamoto/16-00488/QQH4K4F2_S5_L002_R1_001.fastq.gz',
                 u'Shigeru_Miyamoto/16-00489/QQH4K9F2_S6_L002_R1_001.fastq.gz']),
             (u'StephenLavelle', [
                 u'StephenLavelle/16-00982/QQ1H2O1_S1_L001_R1_001.fastq.gz',
                 u'StephenLavelle/16-00983/QQ1H2O2_S2_L001_R1_001.fastq.gz',
                 u'StephenLavelle/16-00984/QQ1H2O3_S3_L001_R1_001.fastq.gz']),
             (u'Undetermined_indices',
              [u'Undetermined_S0_L001_R1_001.fastq.gz',
               u'Undetermined_S0_L002_R1_001.fastq.gz',
               u'Undetermined_S0_L003_R1_001.fastq.gz',
               u'Undetermined_S0_L004_R1_001.fastq.gz',
               u'Undetermined_S0_L005_R1_001.fastq.gz',
               u'Undetermined_S0_L006_R1_001.fastq.gz',
               u'Undetermined_S0_L007_R1_001.fastq.gz',
               u'Undetermined_S0_L008_R1_001.fastq.gz'])])

        self.assertDictEqual(mapping, expected)

        bcl2fastq_output_path = path.join(self.run3_dir,
                                          'Data/Intensities/BaseCalls')

        mapping = get_sample_project_mapping(bcl2fastq_output_path)
        expected = OrderedDict([('',
                                 [u'BUGS-1_CACGTCTA_L001_R1_001.fastq.gz',
                                  u'BUGS-1_CACGTCTA_L001_R2_001.fastq.gz',
                                  u'BUGS-2_AGCTAGTG_L001_R1_001.fastq.gz',
                                  u'BUGS-2_AGCTAGTG_L001_R2_001.fastq.gz',
                                  u'DRUGS-1_ACGTCGTT_L001_R1_001.fastq.gz',
                                  u'DRUGS-1_ACGTCGTT_L001_R2_001.fastq.gz',
                                  u'DRUGS-2_GTCCTGTT_L001_R1_001.fastq.gz',
                                  u'DRUGS-2_GTCCTGTT_L001_R2_001.fastq.gz']),
                                (u'Undetermined_indices',
                                 [u'lane1_Undetermined_L001_R1_001.fastq.gz',
                                  u'lane1_Undetermined_L001_R2_001.fastq.gz'])])

        self.assertDictEqual(mapping, expected)

        # test absolute_paths flag
        mapping = get_sample_project_mapping(bcl2fastq_output_path,
                                             absolute_paths=True)
        for k, v in expected.items():
            expected[k] = [path.join(bcl2fastq_output_path, f) for f in v]

        self.assertDictEqual(mapping, expected)

    # TODO: These times and RTA versions aren't actually consistent
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
