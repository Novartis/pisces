#!/usr/bin/env python
# pylint: disable=R0913, R0914, C0301, C0111, C0103, E0401

import argparse
import os
import sys
import time
import shutil
import json
import socket
import re
import simplesam
import os.path
import pisces.run
import pisces.index
import pisces.qc
import pisces.submit
from signal import signal, SIGINT, SIGTERM
from multiprocessing import Process, cpu_count
from subprocess import Popen, PIPE, call
from collections import OrderedDict, defaultdict
from tqdm import tqdm
from pisces import find_data_directory
from pkg_resources import get_distribution

__version__ = get_distribution("novartis_pisces").version

data_dir = find_data_directory()

class DefaultOrderedDict(OrderedDict):
    def __init__(self, default, items=[]):
        super(DefaultOrderedDict, self).__init__(items)
        self._default = default

    def __missing__(self, key):
        self[key] = value = self._default()
        return value

class ActionRunRscript(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        """ Run an R script """
        from subprocess import call
        import shlex
        call([
            'Rscript',
        ] + [os.path.join(data_dir, 'R/make_expression_matrix.R')] + ['-h'])
        parser.exit()
        
def call_Rscript(args, unknown_args):
    """ Call an R script, passing through arguments from argparse. """
    from tempfile import NamedTemporaryFile
    with NamedTemporaryFile() as config:
        config.write(json.dumps(args.conf).encode())
        config.flush()
        cmd = [
            'Rscript',
            os.path.join(data_dir, 'R/make_expression_matrix.R'), '--config',
            config.name
        ]
        cmd.extend(unknown_args)
        call(cmd)

def default_species_index(conf):
    """ Return a list of the default index builds for each species defined in args.conf """
    defaults = {}
    for species, build in conf.items():
        defaults[species] = build["default"]
    return defaults

def create_parser(args=None):
    # have to hack this: http://bugs.python.org/msg132327

    parser = argparse.ArgumentParser(description="PISCES launch script", add_help=False)
    parser.add_argument('--config', type=str, default=os.path.join(data_dir, 'config.json'), help='default=%(default)s')

    pre, _ = parser.parse_known_args()
    conf = json.loads(open(pre.config).read())
    for species, index_conf in conf.items():
        conf[species]["downloads"] = os.path.expandvars(index_conf["downloads"])
        conf[species]["index_dir"] = os.path.expandvars(index_conf["index_dir"])
    parser.set_defaults(conf=conf)
    
    parser.add_argument('--version', action="version", version=pisces.__version__, help="print PISCES version number")
    parser.add_argument('--debug', action="store_true", help='enable debugging mode')
    parser.add_argument('-h', '--help', action="help")
    

    subparsers = parser.add_subparsers(
        dest="subparser_name",
        help='for help with a subcommand use [subcommand] -h',
        title="valid subcommands",
        metavar='')

    cluster = subparsers.add_parser(
        'submit',
        help=
        "Submit multiple PISCES run jobs to a DRMAA cluster. Extra 'pisces run' options can be passed after the 'metadata' argument."
    )
    cluster.add_argument(
        '--metadata',
        '-m',
        type=argparse.FileType('r'),
        help=
        "metadata.csv file containing (at minimum) SampleID, Fastq1, Fastq2, Directory columns"
    )
    cluster.add_argument(
        '--workdir',
        '-w',
        type=str,
        default=os.getcwd(),
        help=
        "directory where PISCES will create the .pisces folder to contain jobs scripts, logs, and run information"
    )
    cluster.add_argument(
        '--local',
        action="store_true",
        help="run jobs on the local machine")
    cluster.add_argument(
        '--batch',
        action="store_true",
        help="after submitting jobs using DRMAA, exit without monitoring job status")
    cluster.add_argument(
        '--runtime',
        '-rt',
        default=14400,
        type=int,
        help="runtime in seconds for each cluster job")
    cluster.add_argument(
        '--max-memory',
        '-mm',
        default=48,
        type=int,
        choices=range(1, 128),
        help="memory in GB required per job")
    cluster.add_argument(
        '--dry-run',
        '-n',
        action="store_true",
        help="print job commands and then exit")
    cluster.set_defaults(func=pisces.submit._submit_drmaa)

    pipeline = subparsers.add_parser(
        'run',
        help="Run the PISCES pipeline on a single sample")
    required = pipeline.add_argument_group('required arguments')
    optional = pipeline.add_argument_group('optional arguments')
    sra_or_fq = required.add_mutually_exclusive_group(required=True)
    sra_or_fq.add_argument(
        '-fq1',
        type=str,
        nargs='*',
        help="space-separated list of gzipped FASTQ read 1 files")
    sra_or_fq.add_argument(
        '-sra',
        type=str,
        nargs='*',
        help="NCBI sequence read archive accessions in the form SRR#######")
    optional.add_argument(
        '-fq2',
        type=str,
        nargs='*',
        help="space-separated list of gzipped FASTQ read 2 files")

    optional.add_argument(
        '-n',
        '--name',
        type=str,
        help="sample name used in output files. default=auto")
    optional.add_argument(
        '-o',
        '--out',
        type=str,
        help="path to output directory. default=/path/to/$FQ1/PISCES",
    )
    optional.add_argument(
        '-p',
        '--threads',
        type=int,
        default=cpu_count(),
        help="total number of CPU threads to use default=%(default)s")
    optional.add_argument(
        '-t',
        '--sample-type',
        type=str,
        help=
        "type of the library (defined in --config file) default=auto, choices: {0:s}".
        format(str(conf.keys())))
    optional.add_argument(
        '-i',
        '--salmon-indices',
        type=str,
        nargs='*',
        help="salmon indices to use (defined in --config file) defaults={0:s}".
        format(str(default_species_index(conf))))
    optional.add_argument(
        '-s',
        '--salmon-args',
        type=str,
        nargs='*',
        help="extra arguments to pass to salmon (default=%(default)s)")
    optional.add_argument(
        '-l',
        '--libtype',
        type=str,
        choices=('IU', 'ISF', 'ISR', 'A'),
        help=
        "library geometry for Salmon (http://salmon.readthedocs.org/en/latest/salmon.html#what-s-this-libtype) default=auto"
    )
    optional.add_argument(
        '--scratch-dir',
        type=str,
        help="path to scratch directory default='$(--out)'")
    optional.add_argument(
        '--overwrite',
        action='store_true',
        help="overwrite existing files default=False")
    optional.add_argument(
        '--make-bam',
        action='store_true',
        help="make a BAM file for visualization")
    optional.add_argument(
        '--no-salmon', action='store_true', help="do not run salmon")
    optional.add_argument(
        '--no-fastqp',
        action='store_true',
        help="do not generate read-level qc metrics")
    optional.add_argument(
        '--no-vcf', action='store_true', help="do not generate vcf file")
    optional.add_argument(
        '--sra-enc-dir',
        type=str,
        default=None,
        help="path to NCBI SRA project directory for encrypted dbGaP data")
    pipeline.set_defaults(func=pisces.run.run)

    matrix = subparsers.add_parser(
        'summarize-expression',
        help=
        "Compile an expression matrix from multiple individual PISCES runs",
        add_help=False)
    matrix.set_defaults(func=call_Rscript)

    qc = subparsers.add_parser(
        'summarize-qc',
        help="Create QC tables from multiple individual PISCES runs")
    qc.add_argument(
        'directories',
        type=str,
        nargs='*',
        help="directories containing PISCES runs")
    qc.add_argument(
        '--tab',
        type=argparse.FileType('w'),
        default="QC.table.txt",
        help='output file name for tabular data (default=%(default)s)')
    qc.add_argument(
        '--tall',
        type=argparse.FileType('w'),
        default="QC.skinny.txt",
        help=
        'output file name for spotfire tall/skinny format (default=%(default)s)'
    )
    qc.add_argument(
        '-m',
        '--metadata',
        type=argparse.FileType('r'),
        help='metadata csv file')
    qc.add_argument(
        '-f',
        '--fingerprint',
        type=str,
        default=None,
        help=
        'file path at which to create fingerprint probabilities file (default=None)'
    )
    qc.set_defaults(func=pisces.qc.qc_table)

    index = subparsers.add_parser(
        'index', help="Build index files for PISCES using --config file")
    index.add_argument(
        '-p',
        '--threads',
        type=int,
        default=1,
        help="number of CPU threads to use")
    index.add_argument(
        '-i',
        '--indices',
        type=str,
        help=
        "names of indices to build. If none specified, then all defined inidices will be built."
    )
    index.add_argument(
        '--overwrite',
        action='store_true',
        help="overwrite existing index files")
    index.add_argument(
        '--debug',
        action='store_true',
        help="print debugging output to console")
    index.set_defaults(func=pisces.index.build_index)
    
    if len(sys.argv)==1 and not args:
        parser.print_help()
        sys.exit(1)  

    if args:
        args, _ = parser.parse_known_args(args)
    else:
        args, _ = parser.parse_known_args()

    return parser


def clean(*args):
    import sys
    sys.stderr.write("Exiting...")
    sys.exit(0)


def main():
    for sig in (SIGINT, SIGTERM):
        signal(sig, clean)

    parser = create_parser()
    args, unknown_args = parser.parse_known_args()
    args.func(args, unknown_args)


if __name__ == '__main__':
    main()
