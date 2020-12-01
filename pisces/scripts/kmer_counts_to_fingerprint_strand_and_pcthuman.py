usage = '''
#USAGE: python scripts/kmer_counts_to_fingerprint_strand_and_pcthuman.py --kmer-list data/fp_kmers.txt --sample-name SAMPLENAME --out-vcf SAMPLENAME.fastq_fingerprint.vcf --kmer-counts-first-read kmer_counts1.txt --kmer-counts-second-read kmer_counts2.txt 2> SAMPLENAME.mousehuman_and_strand_stats.txt
# Note most so-called "options" are required, specifically: --kmer-list, --kmer-counts-first-read, --kmer-counts-second-read
# (Where kmer_counts are the counts of the kmers in data/fp_kmers.txt, which should be formatted in quads for each SNP like so:

TCAGGATCCATTGAGGACAAA   fwd     1       155735012       T
TTTGTCCTCAATGGATCCTGA   rev     1       155735012       T
TCAGGATCCACTGAGGACAAA   fwd     1       155735012       C
TTTGTCCTCAGTGGATCCTGA   rev     1       155735012       C
TCATCCTCGTTTCTCTTCTGA   fwd     1       179042227       T
TCAGAAGAGAAACGAGGATGA   rev     1       179042227       T
TCATCCTCGTCTCTCTTCTGA   fwd     1       179042227       C
TCAGAAGAGAGACGAGGATGA   rev     1       179042227       C
...

# and where kmer_counts can be generated like so:

FASTQ1="/dlab/NGS/NGS_Data/Novartis/RNA-seq/Xenograft/PTX/CLR6445/CLR6445_ACTTGA_L002_R1_001.fastq.gz";FASTQ2="/dlab/NGS/NGS_Data/Novartis/RNA-seq/Xenograft/PTX/CLR6445/CLR6445_ACTTGA_L002_R2_001.fastq.gz";(zfgrep -h -o -f data/fp_kmers_col1.txt $FASTQ1 > kmer_counts1.txt & zfgrep -h -o -f data/fp_kmers_col1.txt $FASTQ2 > kmer_counts2.txt)
'''
#tested with python 2.7.5
import sys
from math import log10
import optparse

error_rate = 0.03  #error rate probability used in calculation; relatively high
correct_rate = 1 - error_rate
logErrorP = log10(error_rate)
logCorrectP = log10(correct_rate)
logCoinFlip = log10(0.5)


def parseOptions():
    parser = optparse.OptionParser(usage=usage)
    parser.add_option(
        "--fp-kmers", "--kmer-list", help="fingerprinting kmers", default=None)
    parser.add_option(
        "--sample-name",
        "--name",
        "--samp",
        help="sample name; default:TESTSAMP",
        default="TESTSAMP")
    parser.add_option(
        "--out-vcf",
        "--vcf",
        help=
        "fingerprint vcf output file; default: $SAMPLE_NAME.fastq_fingerprint.vcf",
        default=None)
    parser.add_option(
        "--kmer-counts-first-read",
        "--counts1",
        help=
        "file containing kmer counts from first read (should be merged across fastq files if more than one fastq file for R1)",
        default=None)
    parser.add_option(
        "--kmer-counts-second-read",
        "--counts2",
        help=
        "file containing kmer counts from second read (should be merged across fastq files if more than one fastq file for R2)",
        default=None)
    parser.add_option(
        "--stats-out-file",
        "--stats",
        help=
        "file to write mouse/human and strand statistics too; default: stderr",
        default=None)
    return parser


#logFacts = [0]
#cur = 0
#for i in xrange(1,10001):
#    cur += log10(i)
#    logFacts[i] = cur
def count2genprob(ref, alt):
    #ref = max(10000,ref)
    #alt = max(10000,alt)
    #assume 3% error rate; very approximate
    #ref-ref:
    refp = -10 * (logCorrectP * ref + logErrorP * alt
                  )  # + logFacts[ref] + logFacts[alt])
    #alt-alt:
    altp = -10 * (logCorrectP * alt + logErrorP * ref
                  )  # + logFacts[alt] + logFacts[ref])
    #ref-alt:
    hetp = -10 * (logCoinFlip *
                  (ref + alt))  # + logFacts[ref] + logFacts[alt])
    temp = min(refp, altp, hetp)
    gen = '0/0'
    qual = min(altp, hetp) - temp
    if hetp == temp:
        gen = '0/1'
        qual = min(refp, altp) - temp
    elif altp == temp:
        gen = '1/1'
        qual = min(refp, hetp) - temp
    return gen, qual, ((refp - temp), (hetp - temp), (altp - temp))


def main():
    (opts, args) = parseOptions().parse_args()
    if opts.fp_kmers is None or opts.kmer_counts_first_read is None:
        print(
            "--fp-kmers and --kmer-counts-first-read are all required options.  Run with -h for help",
            file=sys.stderr)
        return
    fh = open(opts.fp_kmers)
    loc2kmerlist = {}
    kmer2count = {}
    kmer2filenumcount = {}
    loc2ref = {}
    loc2alt = {}
    for line in fh:
        vals = line.strip().split()
        loc = tuple(vals[2:4])
        if not loc in loc2kmerlist:
            loc2kmerlist[loc] = [vals[0]]
            loc2ref[loc] = vals[-1]
        else:
            loc2kmerlist[loc].append(vals[0])
            loc2alt[loc] = vals[-1]
        kmer2count[vals[0]] = 0
        kmer2filenumcount[vals[0]] = [0] * 2
    fh.close()

    goodlines = 0
    badlines = 0
    foundreads = set()
    fh = open(opts.kmer_counts_first_read)
    for line in fh:
        vals = line.strip().split(":")
        kmer = vals[-1]
        if not kmer in kmer2count:
            badlines += 1
            continue
        if len(vals) > 1:
            foundreads.add(vals[0])
        goodlines += 1
        kmer2count[kmer] += 1
        kmer2filenumcount[kmer][0] += 1
    fh.close()
    if opts.kmer_counts_second_read:
        fh = open(opts.kmer_counts_second_read)
        for line in fh:
            vals = line.strip().split(":")
            kmer = vals[-1]
            if not kmer in kmer2count:
                badlines += 1
                continue
            if len(vals) > 1 and vals[0] in foundreads:
                continue
            goodlines += 1
            kmer2count[kmer] += 1
            kmer2filenumcount[kmer][1] += 1
        fh.close()

    print("Num kmers found:", goodlines, file=sys.stderr)
    print("Num lines not matching kmers:", badlines, file=sys.stderr)

    tothuman = 0
    totmouse = 0
    totfirstreadsense = 0
    totfirstreadantisense = 0
    totsecondreadsense = 0
    totsecondreadantisense = 0
    if opts.stats_out_file is None:
        statsout = sys.stderr
    else:
        statsout = open(opts.stats_out_file, 'w')

    if opts.out_vcf:
        vcfout = open(opts.out_vcf, 'w')
        print(
            '''##fileformat=VCFv4.1
##FILTER=<ID=LowQual,Description="Low quality">
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
##ZgrepGenotyper="analysis_type=ZgrepGenotyper"
##contig=<ID=1,length=249250621>
##contig=<ID=2,length=243199373>
##contig=<ID=3,length=198022430>
##contig=<ID=4,length=191154276>
##contig=<ID=5,length=180915260>
##contig=<ID=6,length=171115067>
##contig=<ID=7,length=159138663>
##contig=<ID=8,length=146364022>
##contig=<ID=9,length=141213431>
##contig=<ID=10,length=135534747>
##contig=<ID=11,length=135006516>
##contig=<ID=12,length=133851895>
##contig=<ID=13,length=115169878>
##contig=<ID=14,length=107349540>
##contig=<ID=15,length=102531392>
##contig=<ID=16,length=90354753>
##contig=<ID=17,length=81195210>
##contig=<ID=18,length=78077248>
##contig=<ID=19,length=59128983>
##contig=<ID=20,length=63025520>
##contig=<ID=21,length=48129895>
##contig=<ID=22,length=51304566>
##contig=<ID=X,length=155270560>
##contig=<ID=Y,length=59373566>
##contig=<ID=MT,length=16569>
##reference=file:///da/onc/sequencing/rna-seq/GRCh37d/GRCh37d.original_genome.normed.fa
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t''' + opts.sample_name,
            file=vcfout)

    for loc in sorted(list(loc2kmerlist.keys()), key=lambda x:9999999999*int(x[0])+int(x[1]) if x[0][0]!='c' else 9999999999*int(x[0][3:])+int(x[1])):
        if len(
                loc2kmerlist[loc]
        ) == 2:  #strand or human/mouse detection instead of SNP detection
            if loc2ref[loc] == 'sense' or loc2alt[loc] == 'sense':
                if loc2ref[loc] == 'sense':
                    sensemer = loc2kmerlist[loc][0]
                    antisensemer = loc2kmerlist[loc][1]
                else:
                    sensemer = loc2kmerlist[loc][1]
                    antisensemer = loc2kmerlist[loc][0]
                totfirstreadsense += kmer2filenumcount[sensemer][0]
                totfirstreadantisense += kmer2filenumcount[antisensemer][0]
                if opts.kmer_counts_second_read:
                    totsecondreadsense += kmer2filenumcount[sensemer][1]
                    totsecondreadantisense += kmer2filenumcount[antisensemer][
                        1]
                #print (loc, "strandedness",kmer2filenumcount[sensemer][0], kmer2filenumcount[sensemer][1], kmer2filenumcount[antisensemer][0], kmer2filenumcount[antisensemer][1],file=sys.stderr)
                continue
            assert loc2ref[loc] == 'human' or loc2ref[loc] == 'mouse'
            if loc2ref[loc] == 'human':
                assert loc2alt[loc] == 'human'
                tothuman += sum([kmer2count[i] for i in loc2kmerlist[loc]])
            if loc2ref[loc] == 'mouse':
                assert loc2alt[loc] == 'mouse'
                totmouse += sum([kmer2count[i] for i in loc2kmerlist[loc]])
            #print (loc, loc2ref[loc], [kmer2count[i] for  i in loc2kmerlist[loc]], file=sys.stderr)
            continue
        refcount = kmer2count[loc2kmerlist[loc]
                              [0]] + kmer2count[loc2kmerlist[loc][1]]
        altcount = kmer2count[loc2kmerlist[loc]
                              [2]] + kmer2count[loc2kmerlist[loc][3]]
        if refcount + altcount == 0:
            continue
        filter = '.'
        if refcount + altcount < 7:
            filter = 'LowQual'
        gen, qual, pl = count2genprob(refcount, altcount)
        if opts.out_vcf:
            print(
                "\t".join([
                    loc[0], loc[1], '.', loc2ref[loc], loc2alt[loc],
                    "%0.2f" % qual, filter, '.', 'GT:AD:DP:GQ:PL', ":".join([
                        gen,
                        str(refcount) + "," + str(altcount),
                        str(refcount + altcount),
                        str(min(int(qual), 99)),
                        ",".join([str(int(i)) for i in pl])
                    ])
                ]),
                file=vcfout)
    print("Human reads:", tothuman, file=statsout)
    print("Mouse reads:", totmouse, file=statsout)
    print(
        "Percent mouse: %0.2f" % (100 * totmouse / float(totmouse + tothuman)),
        file=statsout)
    print("First read sense:", totfirstreadsense, file=statsout)
    print("First read antisense:", totfirstreadantisense, file=statsout)
    print("Second read sense:", totsecondreadsense, file=statsout)
    print("Second read antisense:", totsecondreadantisense, file=statsout)
    print(
        "(>95% indicates first read is sense strand; <5% indicates second read is sense strand; values in between indicate unstranded)",
        file=statsout)
    print(
        "Percent first read sense strand: %0.2f" %
        (100 * (totfirstreadsense + totsecondreadantisense
                ) / float(totfirstreadsense + totfirstreadantisense +
                          totsecondreadsense + totsecondreadantisense)),
        file=statsout)
    if opts.out_vcf:
        vcfout.close()
    statsout.close()


if __name__ == '__main__':
    main()
