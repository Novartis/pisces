from pisces import __version__
from collections import OrderedDict, defaultdict
from pisces import find_data_directory
from subprocess import Popen, PIPE, call
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

def qc_table(args, unknown_args):
    import pandas as pd
    import logging
    import os
    import re
    from tqdm import tqdm
    from itertools import chain
    import simplesam
    import shutil


    logging.basicConfig(
        level=logging.INFO,
        format='\n%(asctime)s %(levelname)-8s %(message)s',
        datefmt='%m-%d %H:%M')

    fingerprint_vcf_pattern = re.compile(r'.*fastq_fingerprint\.vcf$')
    pct_mouse_pattern = re.compile(r'.*pct_human_mouse$')
    counts_pattern = re.compile(r'counts\.txt$')

    sample_directories = DefaultOrderedDict(dict)

    logging.info("Searching for directories.")
    if args.metadata:
        logging.info("Using Directory column from --metadata file.")
        sample_metadata = pd.read_csv(
            args.metadata, dtype={
                'Directory': object,
                'SampleID': object
            })
        directory_search_paths = sample_metadata.Directory
    else:
        directory_search_paths = args.directories
    with tqdm(leave=False, unit='File', unit_scale=True, position=0) as pbar:
        for directory in directory_search_paths:
            logging.info("Seaching %s" % directory)
            for dirpath, dirnames, filenames in os.walk(directory):

                vcf_filename = tuple([
                        _f for _f in (re.search(fingerprint_vcf_pattern, name)
                                      for name in filenames) if _f
                    ])
                pct_human_mouse = tuple([
                        _f for _f in (re.search(pct_mouse_pattern, name)
                                      for name in filenames) if _f
                    ])
                    
                abs_directory = os.path.abspath(dirpath)
                
    
                if any(vcf_filename) and any(pct_human_mouse):
                    if args.metadata:
                        sample_name = sample_metadata[
                            sample_metadata.Directory ==
                            directory].SampleID.tolist()[0]
                    else:
                        sample_name = abs_directory.split('/')[-3]
            
                    try:
                        sample_directories[sample_name]['vcf'] = os.path.join(
                            abs_directory, vcf_filename[0].string)
                        sample_directories[
                            sample_name]['mouse'] = os.path.join(
                                abs_directory, pct_human_mouse[0].string)
                        sample_directories[sample_name][
                            'dirpath'] = os.path.join(dirpath)
                        logging.info("Found {sample_name}".format(**locals()))
                    except IndexError:
                        logging.warn(
                            "Sample {sample_name} is missing files at {abs_directory}.".
                            format(**locals()))
                        logging.debug("{0}".format(
                            sample_directories[sample_name]))
                        del sample_directories[sample_name]
                pbar.update(len(filenames))

    logging.info("Found {0:n} samples.".format(
        len(list(sample_directories.keys()))))

    if args.fingerprint:
        logging.info("Generating fingerprint file.")
        fingerprint_positional_arguments = [
            (sample_name, files['vcf'])
            for sample_name, files in sample_directories.items()
        ]
        cmd = [
            'python',
            os.path.join(data_dir, 'scripts/create_sample_sets_stripped.py'),
            '--pairs-out', args.fingerprint, '--pairs-out-cutoff', '-99999999'
        ]
        cmd.extend([
            ','.join([vcf, name])
            for name, vcf in fingerprint_positional_arguments
        ])
        logging.info("Fingerprinting command: %s", ' '.join(cmd))
        p1 = Popen(cmd, stdout=PIPE, stderr=PIPE)
        returncode = p1.wait()
        if returncode > 0:
            logging.warn("Fingerprinting failed!")
    logging.info("Finished.")