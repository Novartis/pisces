import os
import logging
import pisces
import time
import socket
from multiprocessing import Process
from pisces import find_data_directory, long_substr
from itertools import chain
from pkg_resources import get_distribution

__version__ = get_distribution("novartis_pisces").version

data_dir = find_data_directory()

def run(args, unknown_args):
    """ Main user script to run the PISCES workflow """
    if not args.out:
        output_dir = os.path.join(
            os.path.dirname(os.path.abspath(args.fq1[0])), 'PISCES')
    else:
        output_dir = os.path.abspath(args.out)

    if os.path.exists(output_dir) and not args.overwrite:
        msg = "PISCES output directory %s already exists!" % output_dir
        raise IOError(msg)
    elif os.path.exists(output_dir) and args.overwrite:
        msg = "PISCES output directory %s already exists. Overwriting..." % output_dir
        #shutil.rmtree(output_dir)
        #os.makedirs(output_dir)
    else:
        msg = "Creating PISCES output directory %s." % output_dir
        os.makedirs(output_dir)

    if args.debug:
        level = logging.DEBUG
    else:
        level = logging.INFO

    logging.basicConfig(
        level=level,
        format='%(asctime)s %(name)-8s %(levelname)-8s %(message)s',
        datefmt='%m-%d %H:%M',
        filename=os.path.join(output_dir, 'pisces.log'),
        filemode='w')
    console = logging.StreamHandler()
    console.setLevel(level)
    # set a format which is simpler for console use
    formatter = logging.Formatter('pisces run: %(levelname)-8s %(message)s')
    console.setFormatter(formatter)
    # add the handler to the root logger
    logging.getLogger('').addHandler(console)

    logging.info("PISCES version %s", __version__)
    logging.info("Using config file %s", args.config)
    logging.info(msg)  # Folder creation text before logging init.
    logging.info("PISCES was run with arguments:")
    params = vars(args)
    del params['func']  # remove non-serializable data
    for opt, value in list(params.items()):
        logging.info("%s: %s", opt, str(value))

    # set up the default sample name
    if not args.name:
        if args.fq2:
            sample_name = long_substr(tuple(chain(args.fq1, args.fq2)))
        else:
            sample_name = long_substr(tuple(chain(args.fq1, args.fq1)))
        sample_name = sample_name.split('/')[-1].split('_')[0]
        logging.info("Sample name: %s", sample_name)
        assert len(sample_name) > 0
    else:
        sample_name = args.name

    start_time = time.time()
    logging.info('Running on %s', socket.gethostname())

    if args.scratch_dir:
        scratch_dir = args.scratch_dir
    else:
        scratch_dir = output_dir

    if args.sra:
        fq1, fq2 = ([], [])
        if args.sra_enc_dir:
            download_dir = args.sra_enc_dir
        else:
            download_dir = output_dir
        for accession in args.sra:
            logging.info("Downloading %s from SRA.", accession)
            fq1_d, fq2_d = pisces.download_from_sra(accession, output_dir,
                                                    download_dir)
            fq1.append(fq1_d)
            fq2.append(fq2_d)
            logging.info("Downloaded %s and %s", fq1_d, fq2_d)
        logging.info("SRA downloads finished.")
    else:
        fq1, fq2 = (args.fq1, args.fq2)


    if not args.no_fastqp:
        logging.info('Starting read-level QC.')
        p1 = Process(
            target=pisces.read_level_qc,
            args=(fq1[0], output_dir, sample_name + '_1'))
        p1.start()
        if fq2:
            p2 = Process(
                target=pisces.read_level_qc,
                args=(fq2[0], output_dir, sample_name + '_2'))
            p2.start()
    
    try:
        strand, species = pisces.fingerprint_sample(
            args, fq1, fq2, data_dir, output_dir, sample_name)
    except RuntimeError:
        if not args.sample_type or not args.libtype:
            logging.error("Genotyping failed: please specify --sample-type and --libtype parameters manually.")
            logging.error("Running salmon with --libtype A")
            libtype = 'A'

    if args.libtype:
        libtype = args.libtype
    elif strand > 95 and fq2:
        libtype = 'ISF'
    elif strand < 5 and fq2:
        libtype = 'ISR'
    elif strand < 95 and strand > 5 and fq2:
        libtype = 'IU'
    elif strand > 95 and not fq2:
        libtype = 'SF'
    elif strand < 5 and not fq2:
        libtype = 'SR'
    elif strand < 95 and strand > 5 and not fq2:
        libtype = 'U'

    if args.sample_type:
        sample_type = args.sample_type
    elif species > 99:
        sample_type = 'mouse'
    elif species < 1:
        sample_type = 'human'
    else:
        sample_type = 'xeno'

    if not args.no_salmon:
        logging.info("Starting Salmon quantitation...")
        salmon_indices = args.salmon_indices
        index_dir_base = args.conf[sample_type]['index_dir']
        if not salmon_indices:
            salmon_indices = args.conf[sample_type]['index'].keys()
        for index in salmon_indices:
            try:
                salmon_idx_dir = os.path.join(index_dir_base, sample_type, index, "salmon")                                 
            except KeyError:
                raise RuntimeError(
                    "Salmon index directory {0} not defined in {1}.".format(
                        salmon_idx_dir, args.config))
            else:
                logging.info("Running salmon with index %s", index)
                pisces.run_salmon(fq1, fq2, salmon_idx_dir, args.threads,
                                  libtype, output_dir, data_dir, index,
                                  args.make_bam)

    if not args.no_fastqp:
        p1.join()
        if fq2:
            p2.join()

    logging.info("Analysis complete in %s seconds.",
                 str(time.time() - start_time))
