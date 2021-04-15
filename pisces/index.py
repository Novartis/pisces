import logging
import os
from pisces import find_data_directory
from pkg_resources import get_distribution

__version__ = get_distribution("novartis_pisces").version

def build_index(args, unknown_args):
    from pyfaidx import Fasta
    from twobitreader import TwoBitFile
    import gffutils
    import gffutils.merge_criteria as mc
    import atexit
    import shutil
    from tqdm import tqdm
    from collections import defaultdict
    from pprint import pprint
    from tempfile import NamedTemporaryFile
    from urllib.parse import urlparse
    from urllib.request import urlopen
    from shutil import copyfileobj
    from subprocess import Popen, PIPE, call

    if args.debug:
        log_level = logging.DEBUG
    else:
        log_level = logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s %(levelname)-8s %(message)s',
        datefmt='%m-%d %H:%M')

    config = args.conf
    logging.info("PISCES version %s", __version__)

    for species, datasets in list(config.items()):
        indices = datasets["index"]
        download_dir = datasets["downloads"]
        index_dir_base = datasets["index_dir"]
        if not os.path.exists(download_dir):
            os.makedirs(download_dir)
        for index_name, dataset in list(indices.items()):
            if args.indices and index_name not in args.indices:
                continue
            pprint(dataset, indent=4)
            options = defaultdict(lambda: True)
            options.update(dataset["options"])
            index_dir_path = os.path.join(index_dir_base, species, index_name)
            if os.path.exists(index_dir_path):
                if args.overwrite:
                    logging.warn(
                        "index directory %s already exists! overwriting",
                        index_dir_path)
                    shutil.rmtree(
                        os.path.join(index_dir_path, "transcripts"))
                    shutil.rmtree(os.path.join(index_dir_path, "salmon"))
                else:
                    continue
            os.makedirs(os.path.join(index_dir_path, "transcripts"))
            os.makedirs(os.path.join(index_dir_path, "salmon"))

            transcripts_fasta_file = os.path.join(
                index_dir_path, "transcripts", "transcripts.fa")
                
            with open(transcripts_fasta_file, 'w') as transcripts_fasta:
                ## all of this URI handling should probably use an existing library like
                ## https://github.com/intake/filesystem_spec
                for fasta_loc in dataset["extra_fastas"]:
                    fasta = urlparse(fasta_loc)
                    if fasta.scheme == '':
                        reference = Fasta(fasta.path)
                    elif fasta.scheme.lower() in ('ftp', 'http', 'https'):
                        _fasta_local_path = os.path.join(
                            download_dir, os.path.basename(fasta.path))
                        logging.info("Downloading %s", fasta.geturl())
                        if not os.path.exists(_fasta_local_path):
                            with urlopen(fasta.geturl()) as _fasta:
                                with open(_fasta_local_path,
                                          'wb') as _fasta_local:
                                    copyfileobj(_fasta, _fasta_local)
                                    _fasta_local.flush()
                                    if fasta.path.endswith('gz'):
                                        logging.info("Decompressing %s",
                                                     fasta.geturl())
                                        call(
                                            ' '.join([
                                                'gzip -dc', _fasta_local_path,
                                                '>',
                                                _fasta_local_path.replace(
                                                    ".gz", "")
                                            ]),
                                            shell=True)
                        if _fasta_local_path.endswith("2bit"):
                            logging.info("Converting %s to FASTA format",
                                         fasta.geturl())
                            twobit = TwoBitFile(_fasta_local_path)
                            if not os.path.exists(
                                    _fasta_local_path.replace("2bit", "fa")):
                                with open(
                                        _fasta_local_path.replace(
                                            "2bit", "fa"), 'w') as fasta:
                                    for chrom in twobit.keys():
                                        fasta.write(">%s\n" % chrom)
                                        fasta.write(str(twobit[chrom]) + '\n')
                            reference = Fasta(
                                _fasta_local_path.replace("2bit", "fa"))
                    
                    with open(_fasta_local_path) as extra:
                        logging.info("Adding entries from %s", fasta)
                        for line in extra:
                            transcripts_fasta.write(line)
                            
                for gtf_loc, fasta_loc in zip(dataset["gtfs"],
                                              dataset["fastas"]):
                    gtf = urlparse(gtf_loc)
                    fasta = urlparse(fasta_loc)
                    assembly = os.path.basename(fasta.path)

                    if fasta.scheme == '':
                        reference = Fasta(fasta.path)
                    elif fasta.scheme.lower() in ('ftp', 'http', 'https'):
                        _fasta_local_path = os.path.join(
                            download_dir, os.path.basename(fasta.path))
                        logging.info("Downloading %s", fasta.geturl())
                        if not os.path.exists(_fasta_local_path):
                            with urlopen(fasta.geturl()) as _fasta:
                                with open(_fasta_local_path,
                                          'wb') as _fasta_local:
                                    copyfileobj(_fasta, _fasta_local)
                                    _fasta_local.flush()
                                    if fasta.path.endswith('gz'):
                                        logging.info("Decompressing %s",
                                                     fasta.geturl())
                                        call(
                                            ' '.join([
                                                'gzip -dc', _fasta_local_path,
                                                '>',
                                                _fasta_local_path.replace(
                                                    ".gz", "")
                                            ]),
                                            shell=True)
                        if _fasta_local_path.endswith("2bit"):
                            logging.info("Converting %s to FASTA format",
                                         fasta.geturl())
                            twobit = TwoBitFile(_fasta_local_path)
                            if not os.path.exists(
                                    _fasta_local_path.replace("2bit", "fa")):
                                with open(
                                        _fasta_local_path.replace(
                                            "2bit", "fa"), 'w') as fasta:
                                    for chrom in twobit.keys():
                                        fasta.write(">%s\n" % chrom)
                                        fasta.write(str(twobit[chrom]) + '\n')
                            reference = Fasta(
                                _fasta_local_path.replace("2bit", "fa"))
                        elif fasta.path.endswith('gz'):
                            reference = Fasta(
                                _fasta_local_path.replace(".gz", ""))
                        else:
                            reference = Fasta(_fasta_local_path)

                    if gtf.scheme == '':
                        database_filename = gtf.path + '.db'
                        if os.path.exists(database_filename):
                            logging.info("Loading existing GTF database file.")
                            db = gffutils.FeatureDB(database_filename)
                        else:
                            logging.info(
                                "Creating GTF database file. This will take some time..."
                            )
                            try:
                                db = gffutils.create_db(
                                    gtf.path,
                                    database_filename,
                                    disable_infer_genes=
                                    not options["infer_features"],
                                    disable_infer_transcripts=
                                    not options["infer_features"])
                            except:
                                tmp_db = os.path.join(download_dir, os.path.basename(gtf.path) + '.db')
                                logging.info(
                                    "Unable to create %s, so using %s",
                                    database_filename, tmp_db)
                                if os.path.exists(tmp_db):
                                    logging.info("Loading existing GTF database file.")
                                    db = gffutils.FeatureDB(tmp_db)
                                else:
                                    db = gffutils.create_db(
                                        gtf.path,
                                        tmp_db,
                                        disable_infer_genes=
                                        not options["infer_features"],
                                        disable_infer_transcripts=
                                        not options["infer_features"])
                    elif gtf.scheme.lower() in ('ftp', 'http', 'https'):
                        _gtf_local_path = os.path.join(download_dir,
                                                       os.path.basename(
                                                           gtf.path))
                        logging.info("Downloading %s", gtf.geturl())
                        if not os.path.exists(_gtf_local_path):
                            with urlopen(gtf.geturl()) as _gtf:
                                with open(_gtf_local_path, 'wb') as _gtf_local:
                                    copyfileobj(_gtf, _gtf_local)
                                    _gtf_local.flush()
                                    if gtf.path.endswith('gz'):
                                        logging.info("Decompressing %s",
                                                     gtf.geturl())
                                        call(
                                            ' '.join([
                                                'gzip -dc', _gtf_local_path,
                                                '>',
                                                _gtf_local_path.replace(
                                                    ".gz", "")
                                            ]),
                                            shell=True)
                                        logging.info(
                                            "Creating GTF database file. This will take some time..."
                                        )
                                        db = gffutils.create_db(
                                            _gtf_local_path.replace(".gz", ""),
                                            _gtf_local_path.replace(
                                                ".gz", "") + '.db',
                                            disable_infer_genes=
                                            not options["infer_features"],
                                            disable_infer_transcripts=
                                            not options["infer_features"])
                                    else:
                                        logging.info(
                                            "Creating GTF database file. This will take some time..."
                                        )
                                        db = gffutils.create_db(
                                            _gtf_local_path,
                                            _gtf_local_path + '.db',
                                            disable_infer_genes=
                                            not options["infer_features"],
                                            disable_infer_transcripts=
                                            not options["infer_features"])
                        elif gtf.path.endswith('gz'):
                            logging.info("Loading existing GTF database file.")
                            db = gffutils.FeatureDB(
                                _gtf_local_path.replace(".gz", "") + '.db')
                        else:
                            logging.info("Loading existing GTF database file.")
                            db = gffutils.FeatureDB(_gtf_local_path)

                    # https://github.com/daler/gffutils/issues/56
                    db.execute('ANALYZE features')
                    #if db.count_features_of_type('intron') == 0 and options["unprocessed_transcripts"]:
                        #logging.info("Inferring intronic sequences...")
                        #db.update(db.create_introns())
                    soft_chars = set(('a', 'c', 'g', 't'))

                    if not options["-k"]:
                        k = 31
                    else:
                        k = options["-k"]

                    gene_tx_file = os.path.join(
                        index_dir_path,
                        assembly + "_transcripts_to_genes.txt")
                    gene_annotation = os.path.join(
                        index_dir_path,
                        assembly + "_gene_annotation.txt")
                        
                    def features_to_string(features, fasta_in, masked=True, strand=True):
                        """ 
                        """
                        sequences = []
                        feature_strand = "."
                        for feature in features:
                            feature_strand = feature.strand
                            sequences.append(
                                feature.sequence(
                                    fasta_in, use_strand=strand))
                        # if the transcript is on the reverse strand, reverse order of exons 
                        # before concatenating
                        if feature_strand == "-":
                            sequences = sequences[::-1]
                        seq = ''.join(sequences)
                        mask_count = sum(seq.count(a) for a in soft_chars)
                        if masked:
                            if mask_count > 0:
                                seq = seq.replace(
                                    'a', 'N').replace('t', 'N').replace(
                                        'c', 'N').replace('g', 'N')
                        try:
                            frac_masked = mask_count / len(seq)
                        except ZeroDivisionError:
                            frac_masked = 0
                        return (seq, frac_masked)

                    with open(gene_tx_file, 'w') as gene2tx, open(
                            gene_annotation, 'w') as annotation:
                        logging.info("Making transcripts_to_genes, annotation and FASTA file for %s",
                                     gtf.path)
                        with tqdm(
                                total=db.count_features_of_type('gene'),
                                unit='gene') as pbar:
                            for gene in db.features_of_type('gene'):
                                transcripts = db.children(gene, featuretype='transcript', order_by='start')
                                for transcript in transcripts:
                                    fa_seq, frac_masked = features_to_string(db.children(transcript, 
                                                                                   featuretype='exon', 
                                                                                   order_by='start'), 
                                                                       reference, 
                                                                       masked=options["masked"])
                                    transcripts_fasta.write('>' + transcript['transcript_id'][0] + '\n')
                                    transcripts_fasta.write(fa_seq + '\n')
                                if options["unprocessed_transcripts"]:
                                    exons = db.children(gene, featuretype='exon', order_by='start') 
                                    merged_exons = db.merge(exons, merge_criteria=(mc.seqid, mc.feature_type, mc.overlap_any_inclusive))
                                    introns = db.interfeatures(merged_exons, new_featuretype='intron')                                                     
                                    transcripts_fasta.write('>' + "intronic_" + gene['gene_id'][0] + '\n')
                                    fa_seq, _ = features_to_string(introns, reference, masked=options["masked"])
                                    transcripts_fasta.write(fa_seq + '\n')
                                
                                first_transcript = next(
                                    db.children(
                                        gene,
                                        featuretype='transcript',
                                        order_by='start'))
                                try:
                                    gene_type = first_transcript['gene_type'][0]
                                except KeyError:
                                    gene_type = 'NA'
                                try:
                                    gene_name = first_transcript['gene_name'][0]
                                except KeyError:
                                    gene_name = 'NA'
                                    
                                exons = db.children(gene, featuretype='exon', order_by='start') 
                                merged_exons = db.merge(exons, merge_criteria=(mc.seqid, mc.feature_type, mc.overlap_any_inclusive))
                                
                                annotation.write(
                                    "{gene}\t{type}\t{name}\t{chrom}\t{start}\t{stop}\t{length}\t{frac_masked}\n".
                                    format(
                                        gene=gene['gene_id'][0],
                                        type=gene_type,
                                        name=gene_name,
                                        start=gene.start,
                                        stop=gene.stop,
                                        chrom=gene.chrom,
                                        length=sum(len(exon) for exon in merged_exons),
                                        frac_masked=str(frac_masked)))
                                transcripts = db.children(
                                    gene,
                                    featuretype='transcript',
                                    order_by='start')
                                for transcript in transcripts:
                                    gene2tx.write("{txp}\t{gene}\n".format(
                                        gene=gene['gene_id'][0],
                                        txp=transcript['transcript_id'][0]))
                                pbar.update(1)
                            
                    if options["intergenes"]:
                        for seqid in reference.keys():
                            logging.info("Merging overlapping genes on %s", seqid)
                            merged_genes = db.merge(db.region(seqid=seqid), merge_criteria=(mc.seqid, mc.feature_type, mc.overlap_any_inclusive))
                            with tqdm(unit='intergene features') as pbar:
                                for intergene in db.interfeatures(merged_genes, new_featuretype='intergene'):
                                    transcripts_fasta.write('>' + 'intergene_' + seqid + "_" + str(intergene.start) + ':' + str(intergene.end) + '\n')
                                    fa_seq, _ = features_to_string([intergene], reference, masked=options["masked"], strand=False)
                                    transcripts_fasta.write(fa_seq + '\n')
                                    pbar.update(1)


            # This needs to happen outside of context handler so FASTA file can be closed properly
            logging.info("Making salmon index files for %s",
                         species + '/' + index_name)
            cmd = [
                os.path.join(find_data_directory(), 'redist', 'salmon',
                             'bin', 'salmon'), 'index', '-p',
                str(args.threads), '-k',
                str(k), '-t', transcripts_fasta.name, '-i',
                os.path.join(index_dir_path, "salmon")
            ]
            logging.debug(' '.join(cmd))
            p = Popen(cmd, stderr=PIPE)
            for line in p.stderr:
                line = line.decode()
                if line.endswith('\n'):
                    logging.info(line.rstrip())
                else:
                    logging.info(line)
            logging.info(line)