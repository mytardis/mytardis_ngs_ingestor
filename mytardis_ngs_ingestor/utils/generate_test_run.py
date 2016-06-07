#!/usr/bin/env python

import sys
import os
import shutil
from datetime import datetime
import gzip
from argparse import ArgumentParser
import random
from distutils.spawn import find_executable
from ..illumina.run_info import parse_samplesheet


# TODO: Finish this RunInfo.xml generation function
#       Need to populate the template dict with more variable we
#       don't currently have, including the 'lanes' dictionary, then call
#       write_runinfo_xml somewhere.
#       Add jinja2 to requirements.txt
#
# from jinja2 import Template
#
# runinfo_xml_template = """<?xml version="1.0"?>
# <RunInfo xmlns:xsd="http://www.w3.org/2001/XMLSchema" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" Version="2">
#   <Run Id="{{ run_id }}" Number="{{ run_number }}">
#     <Flowcell>{{ flowcell_id }}</Flowcell>
#     <Instrument>{{ instrument_id }}</Instrument>
#     <Date>{{ date }}</Date>
#     <Reads>
#       {% for lane in lanes %}
#       <Read Number="{{ lane.number }}" NumCycles="{{ lane.cycles }}" IsIndexedRead="{{ lane.is_index }}" />
#       {% endfor %}
#     </Reads>
#     <FlowcellLayout LaneCount="{{ number_of_lanes }}" SurfaceCount="{{ surface_count }}" SwathCount="{{ swath_count }}" TileCount="{{ tile_count }}" />
#     <AlignToPhiX>
#       {% for lane in lanes %}
#       <Lane>{{ lane.number }}</Lane>
#       {% endfor %}
#     </AlignToPhiX>
#   </Run>
# </RunInfo>
# """
#
#
# def write_runinfo_xml(info):
#     runinfo_text = Template(runinfo_xml_template).render(**info)
#     with open(os.path.join(run_path, 'RunInfo.xml'), 'w') as f:
#         f.write(runinfo_text)
#

def add_options(parser):
    parser.add_argument("-s", "--samplesheet",
                        dest="samplesheet_path",
                        type=str,
                        help="The path to Samplesheet.csv",
                        metavar="SAMPLESHEET_PATH")

    parser.add_argument("-x", "--reads-per-sample",
                        dest="reads_per_sample",
                        type=int,
                        default=10,
                        help="The number of reads per sample.",
                        metavar="READS_PER_SAMPLE")

    parser.add_argument("-f", "--fastq-source",
                        dest="fastq_source",
                        type=str,
                        default='',
                        help="The path to a source fastq.gz file that we will"
                             "sample reads from. If specified, --read-length is"
                             "ignored.",
                        metavar="FASTQ_SOURCE")

    parser.add_argument("-o", "--output-path",
                        dest="output_path",
                        type=str,
                        default='.',
                        help="The path where generated files with be written.",
                        metavar="OUTPUT_PATH")

    parser.add_argument("-l", "--read-length",
                        dest="read_length",
                        type=int,
                        default=151,
                        help="The read lengthm overriding any [Read] setting"
                             "in SampleSheet.csv.",
                        metavar="READS_PER_SAMPLE")

    parser.add_argument("-d", "--run-date",
                        dest="run_date",
                        type=str,
                        default='',
                        help="The date in the format YYMMDD",
                        metavar="RUN_DATE")

    parser.add_argument("-i", "--instrument",
                        dest="instrument_id",
                        type=str,
                        default='M04242',
                        help="The instrument ID",
                        metavar="INSTRUMENT_ID")

    parser.add_argument("-c", "--flowcell",
                        dest="flowcell_id",
                        type=str,
                        default='',
                        help="The flowcell ID (FCID)",
                        metavar="FLOWCELL_ID")

    parser.add_argument("-n", "--run-number",
                        dest="run_number",
                        type=int,
                        default=1,
                        help="The run number (eg 14 becomes 0014)",
                        metavar="RUN_NUMBER")

    parser.add_argument("-r", "--random-shuffle",
                        dest="random_shuffle",
                        action='store_true',
                        help="If a source FASTQ file is provided, shuffle the"
                             "nucleotides in place (without shuffling quality"
                             "scores) so we obscure the original source.")

    parser.add_argument("-p", "--paired",
                        dest="paired",
                        action='store_true',
                        help="Generate reads for paired end / mate paired "
                             "(R1/R2).")

    parser.add_argument("--no-create-project-directories",
                        dest="create_project_directories",
                        action='store_false',
                        help="Don't create directory for each project ID. "
                             "Otherwise a project directory in"
                             "(eg SampleProject in SampleSheet.csv) will be "
                             "created, with nested directories for each "
                             "sample ID")

    parser.add_argument("--sample-directory-prefix",
                        dest="sample_directory_prefix",
                        type=str,
                        default='',
                        help="Prefix this to sample directory names. Older "
                             "versions (~1.8.4) of bcl2fastq prefix 'Sample_'"
                             "to the output directory names")

    parser.add_argument("--no-generate-undetermined-indices",
                        dest="generate_undetermined_indices",
                        action='store_false',
                        help="Don't generate Undetermined_indices_*.fastq.gz "
                             "files")

    parser.add_argument("--random-seed",
                        dest="random_seed",
                        type=int,
                        default=-1,
                        help="Set the pseudo-random seed so we can be "
                             "deterministic.",
                        metavar="RANDOM_SEED")


def shuffle_seq(seq):
    l = list(seq)
    random.shuffle(l)
    seq = ''.join(l)
    return seq


def random_seq(length, characters='ATGC'):
    return ''.join([random.choice(characters)
                    for i in range(length)])


def set_missing_key(d, k, v):
    d[k] = d.get(k, v)


def fasta_header(header_type=2, **kwargs):
    # alternative headers in the wild as documented on Wikipedia:
    # https://en.wikipedia.org/wiki/FASTQ_format
    illumina_type_1 = "@{instrument_id}:{lane}:{tile}:{tile_x}:{tile_y}#" \
                      "{index_number}/{pair}"

    illumina_type_2 = "@{instrument_id}:{run_number}:{flowcell_id}:{lane}:" \
                      "{tile}:{tile_x}:{tile_y} {pair}:{filtered}:" \
                      "{control_bits}:{index}"

    illumina_type_3 = "@{instrument_id}:{run_number}:{flowcell_id}:{lane}:" \
                      "{tile}:{tile_x}:{tile_y} {pair}:{filtered}:" \
                      "{control_bits}:{sample_number}"

    fq_header_types = [None, illumina_type_1, illumina_type_2, illumina_type_3]

    set_missing_key(kwargs, 'control_bits', 0)

    return fq_header_types[header_type].format(**kwargs)


def fastq_filename(name_type=1, suffix='.fastq', **kwargs):
    # bcl2fastq 1.8.4 style filenames:
    if kwargs.get('index', None):  # for Undetermined_indices
        template_1 = "{sample_name}_{index}_L{lane:03d}_R{read}_001%s" % suffix
    else:
        template_1 = "{sample_name}_L{lane:03d}_R{read}_001%s" % suffix

    # bcl2fastq 2.x style filenames:
    template_2 = "{sample_name}_{sample_number}_L{lane:03d}_R{read}_001%s" % suffix

    templates = [None, template_1, template_2]

    return templates[name_type].format(**kwargs)


def run_directory_name(instrument_id, run_number, flowcell_id, run_date=None):
    if not run_date:
        run_date = datetime.now().strftime('%y%m%d')

    return "{run_date}_{instrument_id}_{run_number:04d}_{flowcell_id}".format(
        run_date=run_date,
        instrument_id=instrument_id,
        run_number=run_number,
        flowcell_id=flowcell_id)


def write(f, text, to_stdout=True):
    if to_stdout:
        print(text)
    f.write('%s\n' % text)


def run_commandline():
    parser = ArgumentParser(description='Generate FASTQ files that appear'
                                        'to have come from an Illumina '
                                        'instrument run processed by '
                                        'bcl2fastq.')
    add_options(parser)

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    options = parser.parse_args()

    if options.random_seed == -1:
        random.seed()
    else:
        random.seed(options.random_seed)

    number_of_tiles = 5000
    tile_x_max = 200000
    tile_y_max = 200000
    filtered = 'N'
    control_bits = 0
    quality_range = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLM" \
                    "NOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"

    if 'samplesheet_path' not in options:
        sys.stderr.write("Error: no --samplesheet option provided")
        sys.exit(1)

    bgzip_bin = find_executable('bgzip')

    fqsource = None
    if options.fastq_source:
        fqsource = gzip.open(options.fastq_source, 'r')

    samplesheet, chemistry = parse_samplesheet(options.samplesheet_path)

    # make all the SampleSheet keys lowercase
    samplesheet_lowercase = []
    for s in samplesheet:
        r = {k.lower(): s[k] for k in s.keys()}
        samplesheet_lowercase.append(r)
    samplesheet = samplesheet_lowercase

    # if no flowcell id is specified, take the SampleSheet one.
    # if there is no FCID in the SampleSheet, use a default
    if options.flowcell_id:
        flowcell_id = options.flowcell_id
    elif 'FCID' not in samplesheet[0]:
        flowcell_id = '000000000-ANV1L'

    run_id = run_directory_name(options.instrument_id,
                                options.run_number,
                                flowcell_id,
                                run_date=options.run_date)
    run_path = os.path.join(options.output_path, run_id)

    if not os.path.exists(run_path):
        os.makedirs(run_path)

    try:
        shutil.copy(options.samplesheet_path, run_path)
    except shutil.Error:
        # this may occur if the same SampleSheet.csv already exists at this path
        pass

    # TODO: Use [Read] from IEMv4 SampleSheet if --read-length isn't
    #       specified, fallback to default 151 if none specified anywhere
    read_length = options.read_length

    read_pairs = ['1']
    if options.paired:
        read_pairs = ['1', '2']

    # determine the number of lanes based on the highest lane in the SampleSheet
    number_of_lanes = 1
    for s in samplesheet:
        lane = s.get('lane', 0)
        if lane > number_of_lanes:
            number_of_lanes = lane

    # we add stuff to the samplesheet data structure to so that
    # Undetermined_indicies gets generated as if it's a project+samples
    if options.generate_undetermined_indices:
        from collections import OrderedDict

        for pair in read_pairs:
            for lane in range(1, number_of_lanes + 1):
                d = {'lane': lane,
                     'sampleproject': 'Undetermined_indicies',
                     'samplename': 'lane%s_Undetermined' % lane,
                     'sampleid': 'Sample_lane%s' % lane,
                     'pair': pair,
                     'read': pair,
                     'fcid': samplesheet[0].get('fcid', ''),
                     'instrument_id': options.instrument_id,
                     'run_number': options.run_number,
                     }
                samplesheet.append(OrderedDict(d))

    for pair in read_pairs:
        sample_number = 0
        for sample in samplesheet:
            sample_number += 1
            info = {}
            info['pair'] = pair
            info['read'] = pair
            info['project'] = sample.get('sampleproject', '')
            info['sample_number'] = sample_number
            info['sample_name'] = sample.get('samplename', None)
            if info['sample_name'] is None:
                info['sample_name'] = 'Sample-%s' % sample_number

            info['sample_id'] = sample.get('sampleid', '')
            info['lane'] = int(sample.get('lane', 1))
            info['index'] = sample.get('index', '')
            info['flowcell_id'] = sample.get('fcid', '')
            info['instrument_id'] = options.instrument_id
            info['run_number'] = options.run_number

            output_path = os.path.join(run_path, 'Data/Intensities/BaseCalls')

            # progressively build output path, depending on existance of project
            # and sample id values
            if options.create_project_directories and info['project']:
                output_path = os.path.join(output_path, info['project'])
            if options.create_project_directories and info['sample_id']:
                output_path = os.path.join(
                    output_path,
                    "%s%s" % (options.sample_directory_prefix,
                              info['sample_id']))

            if not os.path.exists(output_path):
                os.makedirs(output_path)

            output_filename = fastq_filename(**info)
            fq_outpath = os.path.join(output_path, output_filename)

            with open(fq_outpath, 'w') as fh:
                for n in range(options.reads_per_sample):
                    info['filtered'] = filtered
                    info['control_bits'] = control_bits
                    info['tile'] = random.randrange(1, number_of_tiles)
                    info['tile_x'] = random.randrange(1, tile_x_max)
                    info['tile_y'] = random.randrange(1, tile_y_max)

                    write(fh, fasta_header(**info))

                    if fqsource:
                        # skip header in source
                        try:
                            fqsource.readline()
                        except StopIteration:
                            fqsource = gzip.open(options.fastq_source, 'r')
                        seq = fqsource.readline().strip()
                        if options.random_shuffle:
                            seq = shuffle_seq(seq)
                        plus = fqsource.readline().strip()
                        quality = fqsource.readline().strip()
                    else:
                        seq = random_seq(read_length)
                        plus = '+'
                        # TODO: generating random quality scores is going to
                        #       look silly in FastQC etc - generate some kind of
                        #       plausible quality falloff. Maybe it just makes
                        #       sense to support wrapping ART to generate
                        #       the FASTQ source sequences
                        quality = random_seq(read_length,
                                             characters=quality_range)
                    write(fh, seq)
                    write(fh, plus)
                    write(fh, quality)

            # bgzip them the correct way for FASTQs, falling back to plain gzip
            if bgzip_bin:
                os.system('%s -f %s' % (bgzip_bin, fq_outpath))
            else:
                os.system('gzip -f %s' % fq_outpath)

    if fqsource:
        fqsource.close()

if __name__ == '__main__':
    # Hint - on the commandline, run like:
    #
    #   $ python -m mytardis_ngs_ingestor.utils.generate_test_run
    #
    # (otherwise relative imports will fail)
    run_commandline()
