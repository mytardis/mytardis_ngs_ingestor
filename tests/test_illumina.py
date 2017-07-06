import os
from os import path
import unittest
from datetime import datetime
from collections import OrderedDict
from semantic_version import Version as SemanticVersion
from mytardis_ngs_ingestor.illumina_uploader import \
    get_fastqc_summary_for_project
import mytardis_ngs_ingestor
from mytardis_ngs_ingestor.illumina.run_info import \
    parse_samplesheet, \
    filter_samplesheet_by_project, \
    filter_samplesheet_by_project, \
    rta_complete_parser, get_sample_project_mapping, \
    parse_sample_info_from_filename


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

        self.run5_dir = path.join(
            path.dirname(__file__),
            'test_data/runs/170907_M04242_0005_000000000-ANV1L')

        self.samplesheet_missing_path = path.join(os.path.dirname(__file__),
                                                  self.run5_dir,
                                                  'SampleSheet.csv')

        self.samplesheet_csv_path = path.join(os.path.dirname(__file__),
                                              self.run1_dir,
                                              'SampleSheet.csv')
        self.samplesheet_v4_path = path.join(os.path.dirname(__file__),
                                             self.run2_dir,
                                             'SampleSheet.csv')
        self.samplesheet_single_proj_path = path.join(
            os.path.dirname(__file__),
            'test_data',
            'SampleSheet_single_project.csv')

    def tearDown(self):
        pass

    def test_parse_samplesheet_missing(self):
        samples, chemistry = parse_samplesheet(self.samplesheet_missing_path,
                                               standardize_keys=True,
                                               allow_missing=True)

        self.assertListEqual(samples, [])
        self.assertEqual(chemistry, '')

        with self.assertRaises(IOError) as context:
            samples, chemistry = parse_samplesheet(
                self.samplesheet_missing_path,
                standardize_keys=True,
                allow_missing=False)

    def test_parse_samplesheet_v4(self):
        samples, chemistry = parse_samplesheet(self.samplesheet_v4_path,
                                               standardize_keys=True)
        self.assertEqual(chemistry, 'TruSeq LT')

        expected = [
            {'index': 'CGATGT', 'lane': '1', 'description': '',
             'sampleid': '16-00982', 'sampleplate': '', 'i7indexid': 'A002',
             'samplewell': '', 'sampleproject': 'StephenLavelle',
             'samplename': 'QQ1H2O1'},
            {'index': 'TGACCA', 'lane': '1', 'description': '',
             'sampleid': '16-00983', 'sampleplate': '', 'i7indexid': 'A004',
             'samplewell': '', 'sampleproject': 'StephenLavelle',
             'samplename': 'QQ1H2O2'},
            {'index': 'CAGATC', 'lane': '1', 'description': '',
             'sampleid': '16-00984', 'sampleplate': '', 'i7indexid': 'A007',
             'samplewell': '', 'sampleproject': 'StephenLavelle',
             'samplename': 'QQ1H2O3'},
            {'index': 'AACCAG', 'lane': '2', 'description': '',
             'sampleid': '16-00487', 'sampleplate': '', 'i7indexid': 'A001',
             'samplewell': '', 'sampleproject': 'Shigeru_Miyamoto',
             'samplename': 'QQInputF2'},
            {'index': 'TGGTGA', 'lane': '2', 'description': '',
             'sampleid': '16-00488', 'sampleplate': '', 'i7indexid': 'A002',
             'samplewell': '', 'sampleproject': 'Shigeru_Miyamoto',
             'samplename': 'QQH4K4F2'},
            {'index': 'AGTGAG', 'lane': '2', 'description': '',
             'sampleid': '16-00489', 'sampleplate': '', 'i7indexid': 'A003',
             'samplewell': '', 'sampleproject': 'Shigeru_Miyamoto',
             'samplename': 'QQH4K9F2'},
            {'index': 'AACCAG', 'lane': '3', 'description': '',
             'sampleid': '16-01787', 'sampleplate': '', 'i7indexid': 'A001',
             'samplewell': '', 'sampleproject': 'Phr00t',
             'samplename': 'Q1N'},
            {'index': 'TGGTGA', 'lane': '3', 'description': '',
             'sampleid': '16-01788', 'sampleplate': '', 'i7indexid': 'A002',
             'samplewell': '', 'sampleproject': 'Phr00t',
             'samplename': 'Q1L'},
            {'index': 'AACCAG', 'lane': '4', 'description': '',
             'sampleid': '16-01787', 'sampleplate': '', 'i7indexid': 'A001',
             'samplewell': '', 'sampleproject': 'Phr00t',
             'samplename': 'Q1N'}
        ]

        for expected_sample, sample in zip(expected, samples):
            self.assertDictEqual(sample, expected_sample)

        self.assertEqual(chemistry, 'TruSeq LT')

    def test_parse_samplesheet_csv(self):
        samples, chemistry = parse_samplesheet(self.samplesheet_csv_path,
                                               standardize_keys=True)
        expected = [
            {'control': 'N', 'index': 'AACCAG', 'lane': '1',
             'description': 'May contain nuts',
             'sampleproject': 'GusFring', 'recipe': '',
             'sampleid': '14-06205-OCT4-5', 'fcid': 'H9PJLADXZ',
             'sampleref': 'Hg19', 'operator': 'TW'},
            {'control': 'N', 'index': 'TGGTGA', 'lane': '1', 'description': '',
             'sampleproject': 'GusFring', 'recipe': '',
             'sampleid': '14-06206-OCT4-15', 'fcid': 'H9PJLADXZ',
             'sampleref': 'Hg19', 'operator': 'TW'},
            {'control': 'N', 'index': 'AGTGAG', 'lane': '1', 'description': '',
             'sampleproject': 'GusFring', 'recipe': '',
             'sampleid': '14-06207-ZAX-5', 'fcid': 'H9PJLADXZ',
             'sampleref': 'Hg19', 'operator': 'TW'},
            {'control': 'N', 'index': 'GCACTA', 'lane': '1', 'description': '',
             'sampleproject': 'GusFring', 'recipe': '',
             'sampleid': '14-06208-ZAX-15', 'fcid': 'H9PJLADXZ',
             'sampleref': 'Hg19', 'operator': 'TW'},
            {'control': 'N', 'index': 'TTGGCA', 'lane': '1', 'description': '',
             'sampleproject': 'GusFring', 'recipe': '',
             'sampleid': '14-06200-Input', 'fcid': 'H9PJLADXZ',
             'sampleref': 'Hg19', 'operator': 'TW'},
            {'control': 'N', 'index': 'ACCTCA', 'lane': '2', 'description': '',
             'sampleproject': 'Walter_White', 'recipe': '',
             'sampleid': '14-05655-SW38', 'fcid': 'H9PJLADXZ',
             'sampleref': 'Hg19', 'operator': 'TW'},
            {'control': 'N', 'index': 'AAGAGG', 'lane': '2', 'description': '',
             'sampleproject': 'Walter_White', 'recipe': '',
             'sampleid': '14-05658-SW41', 'fcid': 'H9PJLADXZ',
             'sampleref': 'Hg19', 'operator': 'TW'},
            {'control': 'N', 'index': 'GGAGAA', 'lane': '2', 'description': '',
             'sampleproject': 'Walter_White', 'recipe': '',
             'sampleid': '14-05659-SW42', 'fcid': 'H9PJLADXZ',
             'sampleref': 'Hg19', 'operator': 'TW'},
            {'control': 'N', 'index': 'AGCATG', 'lane': '2', 'description': '',
             'sampleproject': 'Walter_White', 'recipe': '',
             'sampleid': '14-05660-SW43', 'fcid': 'H9PJLADXZ',
             'sampleref': 'Hg19', 'operator': 'TW'},
            {'control': 'N', 'index': 'GAGTCA', 'lane': '2', 'description': '',
             'sampleproject': 'Walter_White', 'recipe': '',
             'sampleid': '14-05661-SW44', 'fcid': 'H9PJLADXZ',
             'sampleref': 'Hg19', 'operator': 'TW'},
            {'control': 'N', 'index': 'CGTAGA', 'lane': '2', 'description': '',
             'sampleproject': 'Walter_White', 'recipe': '',
             'sampleid': '14-06203-SW45', 'fcid': 'H9PJLADXZ',
             'sampleref': 'Hg19', 'operator': 'TW'},
            {'control': 'N', 'index': 'TCAGAG', 'lane': '2', 'description': '',
             'sampleproject': 'Walter_White', 'recipe': '',
             'sampleid': '14-06204-SW46', 'fcid': 'H9PJLADXZ',
             'sampleref': 'Hg19', 'operator': 'TW'}]

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

        project_lines = filter_samplesheet_by_project(
            self.samplesheet_single_proj_path, '', output_ini_headers=True)
        expected = \
"""[Header],,,,,,,\r
IEMFileVersion,4,,,,,,\r
Experiment Name,John Smythe,,,,,,\r
Date,7/09/2015,,,,,,\r
Workflow,GenerateFASTQ,,,,,,\r
Application,FASTQ Only,,,,,,\r
Assay,OTE,,,,,,\r
Description,,,,,,,\r
Chemistry,Default,,,,,,\r
,,,,,,,\r
[Reads],,,,,,,\r
151,,,,,,,\r
151,,,,,,,\r
,,,,,,,\r
[Settings],,,,,,,\r
ReverseComplement,0,,,,,,\r
OnlyGenerateFASTQ,1,,,,,,\r
,,,,,,,\r
[Data],,,,,,,\r
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,Sample_Project,Description\r
1,Sample-04345,BUGS-1,,,A1,CACGTCTA,Cytogenetics,\r
1,Sample-04353,DRUGS-1,,,A2,ACGTCGTT,Cytogenetics,\r
"""

        for expected_line, project_line in \
                zip(expected.splitlines(), project_lines):
            self.assertEqual(project_line, expected_line + '\r\n')



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

    def test_parse_sample_info_from_filename(self):
        fq_info = parse_sample_info_from_filename(
            'DMSO-7_S7_L008_I2_001.fastq.gz')
        self.assertIsNotNone(fq_info)
        self.assertEqual(fq_info.get('sample_name', None), 'DMSO-7')
        self.assertEqual(fq_info.get('sample_number', None), 7)
        self.assertEqual(fq_info.get('lane', None), 8)
        self.assertEqual(fq_info.get('read', None), 2)
        self.assertEqual(fq_info.get('read_type', None), 'I')
        self.assertEqual(fq_info.get('set_number', None), 1)

        fq_info = parse_sample_info_from_filename(
            '16-04333_TGACCA_L007_R3_001.fastq.gz')
        self.assertIsNotNone(fq_info)
        self.assertEqual(fq_info.get('sample_name', None), '16-04333')
        self.assertEqual(fq_info.get('sample_number', None), None)
        self.assertEqual(fq_info.get('lane', None), 7)
        self.assertEqual(fq_info.get('read', None), 3)
        self.assertEqual(fq_info.get('read_type', None), 'R')
        self.assertEqual(fq_info.get('set_number', None), 1)

        fq_info = parse_sample_info_from_filename(
            '//srv/ceph/abcd/ABCD/instrument_runs/nextseq/'
            '170508_NS500295_0103_AHW2JVBGX2_testing/'
            '170508_NS500295_0103_AHW2JVBGX2_testing.bcl2fastq/'
            'SonnyJim/'
            'R10_S10_R1_001.fastq.gz')
        self.assertIsNotNone(fq_info)
        self.assertEqual(fq_info.get('sample_name', None), 'R10')
        self.assertEqual(fq_info.get('sample_number', None), 10)
        self.assertEqual(fq_info.get('lane', None), None)
        self.assertEqual(fq_info.get('read', None), 1)
        self.assertEqual(fq_info.get('read_type', None), 'R')
        self.assertEqual(fq_info.get('set_number', None), 1)

    def test_project_summary_json(self):
        samplesheet_path = self.samplesheet_v4_path
        fastqc_out_path = path.join(path.dirname(__file__), 'test_data/fastqc/')
        samplesheet, chemistry = parse_samplesheet(samplesheet_path,
                                                   standardize_keys=True)
        expected = {u'fastqc_version': u'0.11.3', u'samples': [{u'read': u'R1', u'basic_stats': {'read_length': 51, 'number_of_reads': 2000, 'percent_gc': 51.0}, u'fastqc_report_filename': u'Q1N_S7_L004_R1_001_fastqc.html', u'sample_id': u'Q1N_S7_L004_R1_001', u'illumina_sample_sheet': {}, u'read_type': u'R', u'filename': u'Q1N_S7_L004_R1_001.fastq.gz', u'index': 'AACCAG', u'lane': 4, u'qc_checks': [(u'Basic Statistics', u'PASS'), (u'Per base sequence quality', u'PASS'), (u'Per sequence quality scores', u'WARN'), (u'Per base sequence content', u'PASS'), (u'Per sequence GC content', u'FAIL'), (u'Per base N content', u'PASS'), (u'Sequence Length Distribution', u'PASS'), (u'Sequence Duplication Levels', u'PASS'), (u'Overrepresented sequences', u'PASS'), (u'Adapter Content', u'PASS'), (u'Kmer Content', u'PASS')], u'indexes': {u'index': u'AACCAG'}, u'sample_name': u'Q1N'}]}
        result = get_fastqc_summary_for_project(fastqc_out_path, samplesheet)

        self.assertDictEqual(result, expected)


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
