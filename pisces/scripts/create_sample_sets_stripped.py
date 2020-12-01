import sys
import optparse
import gc
from math import log10, sqrt
from numpy import array
from collections import OrderedDict
badchr = set(['X', 'chrX', 'Y', 'chrY', 'MT',
              'chrM'])  #not used for hetrate score


def parseOptions():
    usage = "Usage: %prog [options] input_vcf_1[,vcfname1] [input_vcf_2[,vcfname2]] [..] [input_vcf_N[,vcfnameN]]"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option(
        "--pairs-out",
        help=
        "Where to write sample1/sample2/p(same) triples, where p(same)>pairs_out_cutoff.  Default: Do not output.  Default pairs_out_cutoff is 0"
    )
    parser.add_option(
        "--pairs-out-cutoff",
        default=0.0,
        type="int",
        help=
        "Value at or below which sample1/sample2/p(same) triples will not be output. Default: 0.0"
    )
    parser.add_option(
        "--hetrate-out",
        help=
        "Where to write hetrate information for each sample.  Default: None")
    parser.add_option(
        "-d",
        "--default-quality",
        default=20,
        type="int",
        help=
        "Default phred-scaled quality score when one is not present in vcf.  Default: 20"
    )
    parser.add_option(
        "--min-quality",
        default=13,
        type="int",
        help=
        "Minimum phred-scaled quality score; genotypes with quality scores less than this will be ignored.  Default: 13"
    )
    parser.add_option(
        "-m",
        "--max-quality",
        default=40,
        type="int",
        help="Maximum phred-scaled quality score to threshold vcf.  Default: 40"
    )
    parser.add_option(
        "--cross-vcf-only",
        default=False,
        action="store_true",
        help=
        "Set if only want to compare one vcf to another, and not within vcf to itself.  (E.g. if comparing to a 'gold dataset'"
    )
    parser.add_option(
        "--first-vcf-only",
        default=False,
        action="store_true",
        help=
        "Set if only want to compare first vcf to other vcfs; useful if comparing against common reference datasets (e.g. CCLE + PTX)"
    )
    parser.add_option(
        "--alt-frequency-file",
        "--alt-allele-frequencies",
        "--aaf",
        default=None,
        help=
        "File containing a mapping from rsIDs to alternate allele frequency. If not specified, program learns aaf on the fly for each snp for each pair of vcf files."
    )
    parser.add_option(
        "-t",
        "--translater",
        default=[],
        action="append",
        help=
        "File containing sample name mapping (e.g. if want to translate from CLD ids to model numbers; will map first column --> second column; may be specified more than once."
    )
    parser.add_option("--verbose", default=False, action="store_true")
    parser.add_option(
        "--ignore-hets",
        "--hom-only",
        default=False,
        action="store_true",
        help=
        "Set if worried about LOH or false hets from SNPchip; only homozygous calls will be considered.  Recommended to use with --default-quality 30 and --force-default-quality."
    )
    parser.add_option(
        "--force-default-quality",
        default=False,
        action="store_true",
        help="Ignore quality scores from vcf and use default quality instead.")
    return parser


def main():
    (opts, args) = parseOptions().parse_args()
    if len(args) == 0:
        print(
            "at least one vcf file is required; see usage with -h",
            file=sys.stderr)
        sys.exit(0)
    if not opts.hetrate_out:
        hetrateout = None
    else:
        hetrateout = open(opts.hetrate_out, 'w')
        print("Sample\tHetRate\tHetRateZ", file=hetrateout)
    rs2aaf = OrderedDict()
    if opts.alt_frequency_file:
        print("Loading alt frequency file", file=sys.stderr)
        aaffile = open(opts.alt_frequency_file)
        for line in aaffile:
            vals = line.strip().split()
            rs2aaf[vals[0]] = float(vals[1])
        aaffile.close()
    sampletranslater = OrderedDict()
    if opts.translater:
        print("Loading sample translater files", file=sys.stderr)
        for file in opts.translater:
            fh = open(file)
            for line in fh:
                if line[0] == '#':
                    continue
                vals = line.strip().split()
                sampletranslater[vals[0]] = vals[1]
            fh.close()
    samplists = []
    samp2genlists = []
    print("Loading vcf files", file=sys.stderr)
    idx = 0
    vcfnames = []
    for vcffilepair in args:
        idx += 1
        if len(vcffilepair.split(',')) > 1:
            vcffile, vcfname = vcffilepair.split(',')
        else:
            vcffile = vcffilepair
            vcfname = str(idx)
        if vcfname in vcfnames:
            print(
                "ERROR:",
                vcfname,
                "is not a unique vcf name. Please assign a unique name to each vcf file, or none at all.  Aborting.",
                file=sys.stderr)
            sys.exit(0)
        else:
            vcfnames.append(vcfname)
        vcf = open(vcffile)
        samplists.append(getInfo(vcf, vcfname))
        samp2genlists.append(
            getGenotypeData(vcf, samplists[-1], opts.default_quality,
                            opts.max_quality, opts.min_quality,
                            opts.force_default_quality))
        vcf.close()
        temp = set()
        for samp in list(samplists[-1].values()):
            if samp in temp:
                print(
                    "ERROR: duplicate samplename found in vcf",
                    vcfname,
                    ":",
                    samp,
                    file=sys.stderr)
                sys.exit(0)
            temp.add(samp)
        print(
            "...loaded",
            len(samp2genlists[-1].samps),
            "samples and",
            len(samp2genlists[-1].snplist),
            "snps from vcf file:",
            vcffile,
            '("' + vcfname + '")',
            file=sys.stderr)
    print("All VCF files loaded", file=sys.stderr)
    results = OrderedDict()
    calculatedhetrates = set()
    for idx in range(len(samplists)):
        if opts.first_vcf_only and idx > 0:
            break
        gens1 = samp2genlists[idx]
        samp1list = gens1.samps
        for idx2 in range(idx, len(samplists)):
            gens2 = samp2genlists[idx2]
            sys.stderr.write('Comparing VCF ' + vcfnames[idx] + ' to VCF ' +
                             vcfnames[idx2] + "...")
            snpinfo, overlappingsnpgenotypes = compileSNPinfo(
                gens1, gens2, rs2aaf, sampletranslater)
            sys.stderr.write(
                ' using ' + str(len(snpinfo.snplist)) + ' overlapping SNPs \n')
            samp2list = gens2.samps
            for samp1 in samp1list:
                if (idx2 == idx and len(samp1list) >= 10) or (
                        not idx in calculatedhetrates and
                    (len(samp1list) + len(samp2list)) > 10):
                    calculatedhetrates.add(idx)
                    if hetrateout:
                        hetrate, hetratez = hetRateZScore(
                            overlappingsnpgenotypes[samp1], snpinfo)
                        if not hetrateout is None:
                            print(samp1, hetrate, hetratez, file=hetrateout)
                if idx2 == idx and opts.cross_vcf_only:
                    continue
                for samp2 in samp2list:
                    if (samp1, samp2) in results or (samp2, samp1) in results:
                        continue
                    if opts.verbose:
                        print(
                            sys.stderr,
                            "Comparing " + samp1 + " from " + vcfnames[idx] +
                            " to " + samp2 + " from " + vcfnames[idx2],
                            end=' ')
                    psame, pdiff, nummatch, nummismatch, numopphom, numhettohom, numhomtohet = probMatch(
                        overlappingsnpgenotypes, samp1, samp2, snpinfo,
                        opts.ignore_hets)
                    if opts.verbose:
                        print(
                            "(psame,pdiff,nummatch,nummismatch,numopphom,numhettohom,numhomtohet):",
                            psame,
                            pdiff,
                            nummatch,
                            nummismatch,
                            numopphom,
                            numhettohom,
                            numhomtohet,
                            file=sys.stderr)
                    score = psame - pdiff
                    results[(samp1, samp2)] = (score, psame, pdiff, nummatch,
                                               nummismatch, numopphom,
                                               numhettohom, numhomtohet)
    if hetrateout:
        hetrateout.close()
    if opts.pairs_out:
        pairsout = open(opts.pairs_out, 'w')
        print(
            "Sample1\tSample2\tp_same\tModel1\tModel2\tpsame\tpdiff\tnummatch\tnummismatch\tnumopphom\tnumhettohom\tnumhomtohet",
            file=pairsout)
        for samplist in samplists:
            for samp1 in sorted(list(samplist.values())):
                for samplist2 in samplists:
                    for samp2 in sorted(list(samplist2.values())):
                        samppair = (samp1, samp2)
                        if samppair in results and results[samppair][0] > opts.pairs_out_cutoff:
                            print(
                                "\t".join([
                                    samp1, samp2,
                                    "%0.1f" % (results[samppair][0]),
                                    sampletranslater[".".join(
                                        samp1.split(".")[:-1])]
                                    if ".".join(samp1.split(".")[:-1]) in
                                    sampletranslater else samp1.split(".")[0],
                                    sampletranslater[".".join(
                                        samp2.split(".")[:-1])]
                                    if ".".join(samp2.split(".")[:-1]) in
                                    sampletranslater else samp2.split(".")[0],
                                    "%0.1f\t%0.1f\t%i\t%i\t%i\t%i\t%i" % tuple(
                                        results[samppair][1:])
                                ]),
                                file=pairsout)
                        samppair = (samp2, samp1)
                        if samppair in results and results[samppair][0] > opts.pairs_out_cutoff:
                            print(
                                "\t".join([
                                    samp1, samp2,
                                    "%0.1f" % (results[samppair][0]),
                                    sampletranslater[".".join(
                                        samp1.split(".")[:-1])]
                                    if ".".join(samp1.split(".")[:-1]) in
                                    sampletranslater else samp1.split(".")[0],
                                    sampletranslater[".".join(
                                        samp2.split(".")[:-1])]
                                    if ".".join(samp2.split(".")[:-1]) in
                                    sampletranslater else samp2.split(".")[0],
                                    "%0.1f\t%0.1f\t%i\t%i\t%i\t%i\t%i" % tuple(
                                        results[samppair][1:])
                                ]),
                                file=pairsout)


class SNPinfo:
    def __init__(self):
        self.snpinfo = OrderedDict()
        self.snplist = []

    def add(self, snp, snpinfo):
        self.snplist.append(snp)
        self.snpinfo[snp] = snpinfo


def compileSNPinfo(gens1, gens2, rs2aaf=None, sampletranslater=None):
    snpinfo = SNPinfo()
    samp1name2id = OrderedDict()
    samp2id2name = OrderedDict()
    if gens1 == gens2:
        multivcf = False
    else:
        multivcf = True
    overlappingsnpgenotypes = OrderedDict()
    for samp in gens1.samps:
        if ".".join(samp.split(".")[:-1]) in sampletranslater:
            samp1name2id[samp] = sampletranslater[".".join(
                samp.split(".")[:-1])]
        else:
            samp1name2id[samp] = samp.split(".")[0]
        overlappingsnpgenotypes[samp] = []
    if multivcf:
        for samp in gens2.samps:
            if ".".join(samp.split(".")[:-1]) in sampletranslater:
                samp2id2name[sampletranslater[".".join(
                    samp.split(".")[:-1])]] = samp
            else:
                samp2id2name[samp.split(".")[0]] = samp
            overlappingsnpgenotypes[samp] = []
    for snp in gens1.snplist:
        if snp not in gens2.snps:
            continue
        gens1overlap = []
        gens2overlap = []
        gens1acount = 0
        gens1bcount = 0
        gens2acount = 0
        gens2bcount = 0
        for samp in gens1.samps:
            temp1 = gens1.get(samp, snp)
            overlappingsnpgenotypes[samp].append(temp1)
            gen1 = temp1.genotype
            if gen1 >= 0:
                gens1acount += 2 - gen1
                gens1bcount += gen1
            if multivcf and samp1name2id[samp] in samp2id2name:
                gen2 = gens2.get(samp2id2name[samp1name2id[samp]],
                                 snp).genotype
                if gen1 >= 0 and gen2 >= 0:
                    gens1overlap.append(gen1)
                    gens2overlap.append(gen2)
        if multivcf:
            for samp in gens2.samps:
                temp2 = gens2.get(samp, snp)
                overlappingsnpgenotypes[samp].append(temp2)
                gen2 = temp2.genotype
                if gen2 >= 0:
                    gens2acount += 2 - gen2
                    gens2bcount += gen2
        gens1count = float(gens1acount + gens1bcount)
        if gens1count == 0:
            for samp in list(overlappingsnpgenotypes.values()):
                samp.pop()
            continue
        numsamps = len(gens1.samps)
        if multivcf:
            numsamps += len(gens2.samps)
            gens2count = float(gens2acount + gens2bcount)
            if gens2count == 0:
                for samp in list(overlappingsnpgenotypes.values()):
                    samp.pop()
                continue
            baf = (gens1bcount + gens2bcount) / (gens1count + gens2count)
            if numsamps < 10:
                baf = (10 + gens1bcount + gens2bcount) / (
                    20 + gens1count + gens2count)
        else:
            baf = gens1bcount / gens1count
            if numsamps < 10:
                baf = (10 + gens1bcount) / (20 + gens1count)
        if snp[2] in rs2aaf:
            baf = rs2aaf[snp[2]]
        elif not rs2aaf == OrderedDict():
            print(
                'Warning: rs2aaf specified, but snpid',
                snp[2],
                'not found in file',
                file=sys.stderr)
        baf = min(0.99, max(baf, 0.01))
        snpinfo.add(snp, (baf, log10(
            (1 - baf)**2), log10(2 * (1 - baf) * baf), log10(baf**2),
                          (1 - baf)**2, 2 * (1 - baf) * baf, baf**2))
    return snpinfo, overlappingsnpgenotypes


phred2error = [10**(-Q / 10.0) for Q in range(100)]
phred2logprobcorrect = [-10] + [
    log10(1.0 - phred2error[Q]) for Q in range(1, 100)
]


def probMatch(gens, samp1, samp2, snpinfos, ignorehets):
    #TODO / Note: This code may break down if error rates aren't generally small!
    psame = 0
    pdiff = 0
    nummatch = 0
    nummismatch = 0
    numopphom = 0
    numhettohom = 0
    numhomtohet = 0
    idx = 0
    gens1 = gens[samp1]
    gens2 = gens[samp2]
    numsnps = len(gens1)
    while idx < numsnps:
        gen1 = gens1[idx].genotype
        qual1 = gens1[idx].quality
        gen2 = gens2[idx].genotype
        qual2 = gens2[idx].quality
        if gen1 < 0 or gen2 < 0 or (ignorehets and (gen1 == 1 or gen2 == 1)):
            idx += 1
            continue
        assert gen1 <= 2 and gen2 <= 2
        snpinfo = snpinfos.snpinfo[snpinfos.snplist[idx]]
        if (gen1 == gen2):
            #psame ~ p(gen) * p(neither is an error)
            #psame += log10(genprob(snp,gen1[0],snpinfo)) + log10(1-10**(-gen1[1]/10.0)) + log10(1-10**(-gen2[1]/10.0))
            psame += snpinfo[1 +
                             gen1] + phred2logprobcorrect[qual1] + phred2logprobcorrect[qual2]
            nummatch += 1
        else:
            #psame ~ p(gen2)*p(gen1 error) + p(gen1)*p(gen2 error)
            #psame += log10( genprob(snp,gen2[0],snpinfo)*(10**(-gen1[1]/10.0)) +
            #                genprob(snp,gen1[0],snpinfo)*(10**(-gen2[1]/10.0)) )
            psame += log10(snpinfo[4 + gen2] * phred2error[qual1] +
                           snpinfo[4 + gen1] * phred2error[qual2])
            nummismatch += 1
            if gen1 == 1:
                numhettohom += 1
            elif gen2 == 1:
                numhomtohet += 1
            else:
                numopphom += 1
                #print >>sys.stderr, snpinfos.snplist[idx]
        #pdiff = p(gen1) * p(gen2)
        #pdiff += log10(genprob(snp,gen1[0],snpinfo)) + log10(genprob(snp,gen2[0],snpinfo))
        pdiff += snpinfo[1 + gen1] + snpinfo[1 + gen2]
        idx += 1
    return psame, pdiff, nummatch, nummismatch, numopphom, numhettohom, numhomtohet


#def genprob(snp, gen, snpinfo):
#    info = snpinfo[snp]
#    q = info[4]
#    p = 1-q
#    if gen == 0:
#        return p**2
#    if gen == 1:
#        return 2*p*q
#    if gen == 2:
#        return q**2
#    print >>sys.stderr, "Error: unknown genotype "  + str(gen)
#    exit(0)


def hetRateZScore(gens, snpinfos):
    m = 0  #expected number of hets
    v = 0  #expected variance on number of hets
    obs = 0  #number of observed hets
    tot = 0  #number of genotypes
    m2 = 0  #expected number of hets if we are not observing any homref
    v2 = 0  #expected variance on number of hets if we are not observing any homref
    numsnps = len(snpinfos.snplist)
    for idx in range(numsnps):
        chr = snpinfos.snplist[idx][0]
        if chr in badchr:
            continue
        info = snpinfos.snpinfo[snpinfos.snplist[idx]]
        gen = gens[idx].genotype
        if gen < 0 or gens[idx].quality < 20 or info[5] < 0.1:  #if nocall, or low quality, or not a 10% chance of observing het, skip--too much bad data in this range
            continue
        if gen == 1:
            obs += 1
        #q = info[0]
        #p = 1-q
        #phet = 2*p*q
        phet = info[5]
        m += phet
        v += phet * (1 - phet)
        tot += 1
    if tot < 10:
        return None, None
    return float(obs) / tot, (obs - m) / sqrt(v)


def getInfo(vcf, vcfname=None):
    col2samp = OrderedDict()
    samps = set()
    while True:
        line = vcf.readline().strip()
        if line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            vals = line.strip().split()
            for i in range(9, len(vals)):
                if vals[i] in samps:
                    suffixnum = 2
                    suffix = "." + str(suffixnum)
                    while vals[i] + suffix in samps:
                        suffixnum += 1
                        suffix = "." + str(suffixnum)
                    print(
                        "Warning: renaming duplicate sample (within single vcf)",
                        vals[i],
                        vals[i] + suffix,
                        file=sys.stderr)
                    vals[i] = vals[i] + suffix
                samps.add(vals[i])
                if vcfname:
                    col2samp[i] = vals[i] + "." + vcfname
                else:
                    col2samp[i] = vals[i]
            break
    return col2samp


class Genotypes:
    def __init__(self, samps):
        self.samps = list(samps.values())
        self.samps.sort()
        self.col2samp = samps
        self.genotypes = []
        self.samp2col = OrderedDict()
        self.snps = set()
        self.snplist = list(self.snps)
        self.snp2col = OrderedDict()
        for key, value in list(samps.items()):
            self.samp2col[value] = key - 9
        self.snp2col = OrderedDict()

    def add(self, snp, genotypes):
        self.genotypes.append(genotypes)
        self.snps.add(snp)
        self.snp2col[snp] = len(self.genotypes) - 1
        self.snplist = list(self.snps)

    def get(self, samp, snp):
        return self.genotypes[self.snp2col[snp]][self.samp2col[samp]]

    def getlist(self, samp, snplist):
        sampcol = self.samp2col[samp]
        return [self.genotypes[self.snp2col[snp]][sampcol] for snp in snplist]


class Genotype:
    def __init__(self, genotype, quality):
        self.genotype = genotype
        self.quality = quality


def getGenotypeData(vcf, samps, defaultQuality, maxQuality, minQuality,
                    forceDefaultQuality):
    gendata = Genotypes(samps)
    gc.disable()
    for line in vcf:
        vals = line.strip().split()
        if vals[6] != "." and vals[6] != "PASS":
            continue
        snp = tuple(vals[:2] + vals[3:5])
        format = vals[8].split(":")
        GT = -1
        GQ = -1
        for idx, formatstring in enumerate(format):
            if (formatstring == "GT"):
                GT = idx
            elif (formatstring == "GQ"):
                GQ = idx
        assert GT >= 0
        hasGQ = GQ >= 0
        genotypes = [None] * len(samps)
        for i in range(9, len(vals)):
            temp = vals[i].split(":")
            if (hasGQ and len(temp) > GQ):
                gq = int(float(temp[GQ]))
                if gq < minQuality:
                    genotypes[i - 9] = Genotype(gen2val('./.'), gq)
                elif forceDefaultQuality:
                    genotypes[i - 9] = Genotype(
                        gen2val(temp[GT]), defaultQuality)
                else:
                    genotypes[i - 9] = Genotype(
                        gen2val(temp[GT]), min(gq, maxQuality))
            else:
                genotypes[i - 9] = Genotype(gen2val(temp[GT]), defaultQuality)
        gendata.add(snp, genotypes)
    gc.enable()
    return gendata


def gen2val(gen):
    #if(gen=="./."):
    #    return -1
    if (gen == "0/0" or gen == "0|0"):
        return 0
    if (gen == "0/1" or gen == "1/0" or gen == "0|1" or gen == "1|0"):
        return 1
    if (gen == "1/1" or gen == "1|1"):
        return 2
    return -1
    #print >>sys.stderr, "Unrecognizable genotype: " + gen
    #exit(0)


if __name__ == "__main__":
    sys.exit(main())
