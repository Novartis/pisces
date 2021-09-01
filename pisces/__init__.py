# pylint: disable=R0913, R0914, C0301, C0111, C0103

import logging
import os
import re
import sys
import json
import time
import pwd
import atexit
import string
import random
import shutil
import platform
import tarfile
from pprint import pformat
from subprocess import Popen, PIPE, call
from multiprocessing import Process
from tempfile import NamedTemporaryFile, mkdtemp
from functools import partial
from pkg_resources import get_distribution

__version__ = get_distribution("novartis_pisces").version

unique_id = ''.join(random.choice(string.digits) for _ in range(10))    

def find_data_directory():
    """ Returns the path to the module directory """
    return os.path.dirname(os.path.abspath(__file__))

def install_salmon():
    """Install a version of salmon from bioconda"""
    import glob
    import requests
    from io import BytesIO
    from urllib.request import urlopen
    from shutil import rmtree
    from subprocess import call
    
    redist = os.path.join(find_data_directory(), 'redist')
    rmtree(os.path.join(redist, "salmon"), ignore_errors=True)
   
    if platform.system() == "Linux":
        salmon_url = "https://anaconda.org/bioconda/salmon/1.3.0/download/linux-64/salmon-1.3.0-hf69c8f4_0.tar.bz2"
    elif platform.system() == "Darwin":
        salmon_url = "https://anaconda.org/bioconda/salmon/1.3.0/download/osx-64/salmon-1.3.0-hb70dc8d_0.tar.bz2"

    print("Installing salmon")
    with requests.get(salmon_url, allow_redirects=True) as response:
        tar_file = BytesIO(response.content)
        tar_file.seek(0)
        with tarfile.open(fileobj=tar_file) as tar:
            tar.extractall(path=os.path.join(redist, "salmon"))
            
def install_r_dependencies():
    """Install R dependencies"""
    cmd = [
        'Rscript',
        os.path.join(find_data_directory(), 'R/set_up_dependencies.R')
    ]
    call(cmd)
    
def install_dependencies():
    install_salmon()
    install_r_dependencies()

def sra_valid_accession(accession):
    """ Test whether a string is an SRA accession """
    if accession.startswith('SRR') and len(accession) == 10:
        return True
    return False
    
def long_substr(data):
    """ http://stackoverflow.com/questions/2892931/longest-common-substring-from-more-than-two-strings-python """
    substr = ''
    if len(data) > 1 and len(data[0]) > 0:
        for i in range(len(data[0])):
            for j in range(len(data[0]) - i + 1):
                if j > len(substr) and all(data[0][i:i + j] in x
                                           for x in data):
                    substr = data[0][i:i + j]
    return substr

def download_from_sra(accession, output_dir, enc_dir):
    """ Dump FASTQ files to output_dir location """
    fastq1_name = os.path.join(output_dir, accession + '_1.fastq')
    fastq2_name = os.path.join(output_dir, accession + '_2.fastq')

    if os.path.exists(fastq1_name) or os.path.exists(fastq2_name):
        logging.info("SRA download for %s already exists.", accession)
    else:
        os.chdir(enc_dir)
        cmd = [
            'fastq-dump', '--split-spot', '--split-files', '-F',
            '--defline-seq', '@$ac.$si.$sg/$ri', '--defline-qual', '+', '-O',
            os.path.join(output_dir), accession
        ]
        logging.info('fastq-dump command: %s', ' '.join(cmd))
        p1 = Popen(cmd)
        r1 = p1.wait()
        if r1:
            logging.warn("SRA download failed.")
            sys.exit(1)
    if os.path.exists(fastq1_name) and os.path.exists(fastq2_name):
        return (fastq1_name, fastq2_name)
    elif os.path.exists(fastq1_name) and not os.path.exists(fastq2_name):
        return (fastq1_name, None)

def read_level_qc(fastq, output_dir, sample_name):
    from fastqp.cli import run
    # run read-level QC step
    arguments = {
        'quiet': False,
        'binsize': False,
        'nreads': 2000000,
        'base_probs': '0.25,0.25,0.25,0.25,0.1',
        'kmer': 5,
        'type': None,
        'leftlimit': 1,
        'rightlimit': -1,
        'median-qual': 30,
        'aligned_only': False,
        'unaligned_only': False,
        'count_duplicates': False,
        'median_qual': 20
    }

    arguments.update(
        [('input', fastq), ('output',
                            os.path.join(output_dir, sample_name + '_fastqp')),
         ('text', os.path.join(output_dir, sample_name + '_fastqp.txt')),
         ('name', sample_name)])

    logging.info('Starting fastqp analysis for %s', sample_name)
    run(arguments)

def fingerprint_sample(args, fastq_1, fastq_2, data_dir, output_dir,
                       sample_name):
    # run zgrep fingerprinting
    # TODO: actual fingerprinting using create_sample_sets_stripped.py
    logging.info("Started genotyping.")
    start = time.time()
    kmers1_outname = os.path.join(output_dir,
                                  sample_name + '.fastq1_kmers.txt')
    kmers2_outname = os.path.join(output_dir,
                                  sample_name + '.fastq2_kmers.txt')
    if fastq_2:
        with open(kmers1_outname, 'w') as kmers1_out, open(
                kmers2_outname, 'w') as kmers2_out:
            p1 = Popen(
                [
                    'fgrep', '-h', '-o', '-f',
                    os.path.join(data_dir, 'data/fp_kmers_col1.txt'),
                    ' '.join(fastq_1)
                ],
                stdout=kmers1_out,
                stderr=PIPE)
            p2 = Popen(
                [
                    'fgrep', '-h', '-o', '-f',
                    os.path.join(data_dir, 'data/fp_kmers_col1.txt'),
                    ' '.join(fastq_2)
                ],
                stdout=kmers2_out,
                stderr=PIPE)
            p1.wait()
            p2.wait()
    else:
        with open(kmers1_outname, 'w') as kmers1_out:
            p1 = Popen(
                [
                    'fgrep', '-h', '-o', '-f',
                    os.path.join(data_dir, 'data/fp_kmers_col1.txt'),
                    ' '.join(fastq_1)
                ],
                stdout=kmers1_out,
                stderr=PIPE)
            p1.wait()
    logging.info("Kmer counts completed.")
    if not args.no_vcf:
        logging.info("Making VCF")
    stats_file = os.path.join(output_dir, sample_name + '.pct_human_mouse')
    cmd = [
        'python',
        os.path.join(
            data_dir,
            'scripts/kmer_counts_to_fingerprint_strand_and_pcthuman.py'),
        '--fp-kmers',
        os.path.join(data_dir, 'data/fp_kmers.txt'), '--sample-name',
        sample_name, '--counts1', kmers1_outname, '--stats', stats_file
    ]
    if fastq_2:
        cmd.append('--counts2')
        cmd.append(kmers2_outname)
    if not args.no_vcf:
        cmd.append('--out-vcf')
        cmd.append(
            os.path.join(output_dir, sample_name + '.fastq_fingerprint.vcf'))
    logging.debug("Fingerprint command: %s", ' '.join(cmd))
    p3 = Popen(cmd, stderr=PIPE)
    _, err3 = p3.communicate()
    strand = None
    species = None
    with open(stats_file) as stats:
        pattern = re.compile(r'Percent ([a-z ]+): ([0-9.]+)')
        for line in stats:
            match = re.match(pattern, line.rstrip())
            if match:
                if match.group(1).split()[0] == 'mouse':
                    species = float(match.group(2))
                elif match.group(1).split()[0] == 'first':
                    strand = float(match.group(2))
    logging.info("Percent mouse: %s", str(species))
    logging.info("Percent first read sense strand: %s", str(strand))
    logging.debug(err3)
    if strand is not None and species is not None:
        logging.info("Finished genotyping.")
        return (strand, species)
    else:
        logging.warn("Genotyping failed!")
        raise RuntimeError("Genotyping failed.")
    logging.info("Sample identification finished in %s seconds.",
                 time.time() - start)


def run_salmon(fastq_1, fastq_2, index_dir_path, threads, libtype, output_dir,
               data_dir, index_name, make_bam):
    import shutil

    if os.path.exists(os.path.join(output_dir, 'salmon', index_name)):
        logging.info(
            "Salmon output directory exists. Removing existing output.")
        shutil.rmtree(os.path.join(output_dir, 'salmon', index_name))
    os.makedirs(os.path.join(output_dir, 'salmon', index_name))

    with open(
            os.path.join(
                os.path.join(output_dir, 'salmon', index_name),
                "pisces_params.json"), 'w') as pisces_params:
        pisces_params.write(
            json.dumps({
                "index_dir": index_dir_path,
                "index_name": index_name
            }))
    p = Popen(
        format_salmon_command(
            libtype,
            threads,
            # This is a KV pair from config.json
            index_dir_path,
            os.path.join(output_dir, 'salmon', index_name),
            fastq_1,
            fastq_2,
            data_dir,
            make_bam),
        shell=True,
        stderr=PIPE)
    for line in p.stderr:
        line = line.decode()
        if line.endswith('\n'):
            logging.info(line.rstrip())
        else:
            logging.debug(line)


def format_salmon_command(libtype, threads, index, output_dir, read1, read2,
                          data_dir, make_bam):
    if read2:
        cmd = [
            os.path.join(data_dir, 'redist', 'salmon',
                         'bin', 'salmon'), '--no-version-check', 'quant', '-q',
            '--index', index, '--libType', libtype, '--mates1', ' '.join(read1),
            '--mates2', ' '.join(read2), '--output', output_dir, '--threads',
            str(threads), '--seqBias', '--gcBias', '--validateMappings', '--dumpEq'
        ]
    elif not read2:
        cmd = [
            os.path.join(data_dir, 'redist', 'salmon',
                         'bin', 'salmon'), '--no-version-check', 'quant', '-q',
            '--index', index, '--libType', libtype, '-r', ' '.join(read1), '--output',
            output_dir, '--threads',
            str(threads), '--seqBias', '--gcBias', '--validateMappings', '--dumpEq'
        ]
    if make_bam:
        import shlex
        cmd.extend(
            shlex.split(
                "--writeMappings | {0} view -Sb - | {0} sort -T {1} -o -  > {2}".
                format("samtools",
                    os.path.join(output_dir, "sort.tmp"),
                    os.path.join(output_dir, "mapped.bam"))))
    logging.debug("Salmon command line: " + ' '.join(cmd))
    return ' '.join(cmd)
