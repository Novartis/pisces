=============================================================================================
PISCES: a package for rapid quantitation and quality control of large scale mRNA-seq datasets
=============================================================================================
-------------------------------------------------------------------------------------------------------------------
Matthew D. Shirley\ :sup:`1`, Viveksagar K. Radhakrishna\ :sup:`1`, Javad Golji\ :sup:`1`, Joshua M. Korn\ :sup:`1`
-------------------------------------------------------------------------------------------------------------------

Novartis Institutes for Biomedical Research\ :sup:`1`

Contact: matt_d.shirley@novartis.com, joshua.korn@novartis.com

Abstract
========

PISCES eases processing of large mRNA-seq experiments by encouraging capture of
metadata using simple textual file formats, processing samples on either a single machine or in parallel
on a high performance computing cluster (HPC), validating sample identity using
genetic fingerprinting, and summarizing all outputs in analysis-ready data
matrices. PISCES consists of two modules: 1) compute cluster-aware
analysis of individual mRNA-seq libraries including species detection, SNP
genotyping, library geometry detection, and quantitation using salmon, and 2)
gene-level transcript aggregation, transcriptional and read-based QC, TMM normalization and differential expression
analysis of multiple libraries to produce data ready for visualization and
further analysis.

PISCES is implemented as a python3 package and is bundled with all necessary
dependencies to enable reproducible analysis and easy deployment. JSON
configuration files are used to build and identify transcriptome indices, and
CSV files are used to supply sample metadata and to define comparison groups for
differential expression analysis using DEseq2. PISCES builds on many existing
open-source tools, and releases of PISCES are available on GitHub_ or the python
package index (PyPI_).

.. _GitHub: https://github.com/Novartis/pisces
.. _PyPI: https://pypi.org/project/novartis-pisces


Introduction
============

Since the first description of RNA-seq :cite:`mortazavi2008mapping`, methods for
RNAseq quantification have rapidly increased in sensitivity and decreased in
required processing time. Recent improvements in speed have been achieved by
removing the step of full read alignment to a reference
:cite:`Patro2017,Bray2016,Patro2014`. While these “alignment-free” (also known
as “pseudoalignment,” “quasi-mapping,” or “lightweight alignment”) methods have
achieved real speed gains, this comes at a loss of compatibility with existing
RNAseq workflows, including important integrity checks. Alignment-free quality
control (QC) and sample identification methods for RNAseq libraries are not
readily available, and existing packages :cite:`picard` for QC are incompatible
with the data formats produced by these new quantification methods. Decreased
sequencing cost have also led to higher throughput of sequencing, with an
associated desire for faster quantification and turnaround of primary analysis.
Typical RNAseq experiments may now be composed of hundreds to thousands of
individual samples, each with descriptive variables for downstream analysis,
which increases the complexity of data and metadata management. Increased number
of samples presents more opportunities for sample mixups in the lab. Complex
designs or cross-experiment meta-analyses further present opportunities for in
silico label swapping. Tracking of these crucial metadata are usually left to
the sequencing facility or individual analyst. With no widely agreed upon
standards for metadata exchange, information about RNAseq libraries may be lost.

With these types of experiments in mind, we built PISCES with an alignment-free
human SNP fingerprinting method for checking sample identity, efficient FASTQ
QC, and novel QC of the transcriptomic and genomic compartments without the need
for full genomic alignment. PISCES is driven by a well-defined CSV metadata
format (:numref:`sra_table`), which limits opportunities for label swaps and encourages analysts to
retain a minimal set of metadata necessary to reproduce an analysis.

In addition to streamlining primary analysis (transcript and gene
quantification) at the single sample level, PISCES provides methods for
executing experiment-level analyses. One common secondary analysis method for
RNAseq libraries is differential gene expression analysis (DGE). PISCES
implements DGE using a wrapper for DEseq2 :cite:`Love2014`, 
with contrast groups defined using descriptive variables referenced in the CSV
metadata file. PISCES creates summarized data matrices of transcript-level
counts and TPM, as well as gene-level counts, TPM, and log\ :subscript:`2`
fold-change calculated using the median of a user-specified reference group for
each condition. Trimmed mean of M-values normalization :cite:`Robinson2010` as
implemented in edgeR is used to adjust TPM abundances of protein coding genes to
account for differences in RNAseq library transcript composition, e.g. a varying
amount of non-coding transcripts due to pre-mRNA contamination or incomplete
rRNA depletion.

Finally, with increasing numbers of RNAseq libraries comes an increased
computational burden. PISCES communicates directly with modern compute clusters
using the distributed resource management application API (DRMAA) to efficiently
submit and monitor jobs processing individual RNAseq libraries.


Implementation of the PISCES package
====================================

PISCES is provided as a Python 3 package, including all tools and dependencies,
which users can easily install
using common packaging tools such as ``pip`` on modern Linux and MacOS machines.

Versioned releases are available from both
Github_ and  PyPI_. Portions of PISCES are implemented
in the R statistical computing language, and these dependencies, as well as
other binary dependencies, are also bundled 
and installed automatically using renv :cite:`Ushey2020`. Workflows can be run
either locally, on a single
machine, or  on a compute cluster supporting the DRMAA interface. Individual
RNAseq libraries may be analyzed without associated metadata by passing FASTQ
files directly into PISCES, or entire experiments may be defined using a
CSV format metadata file containing, at a minimum: sample names (*SampleID*),
FASTQ file locations (*Fastq1 and optionally Fastq2*) or NCBI SRA run accessions
(*SRA*), and output directories for each sample (*Directory*). The
metadata file may also contain any descriptive information associated with the
biological samples such as treatment, batch, or physical sample QC metrics.
This descriptive information can be used by PISCES to identify sample groupings
and reference samples for calculating normalized log2 fold changes, and to
identify groups of samples for generating QC and exploratory figures.
Differential expression is performed using DEseq2, using contrasts defined in a
separate CSV format (:numref:`contrasts_table`) that describes covariates of
interest referenced in the metadata CSV file. Transcriptomic QC metrics are
calculated from the ratios of unprocessed intronic transcripts, mature processed
transcripts, and intergenic regions. Several library QC metrics are inferred
using only kmer counting directly from FASTQ files, including species detection,
strand detection, and human SNP fingerprint for sample identification.

The PISCES workflow
-------------------

.. _workflow_graph:
.. graphviz::
  :caption: Overview of PISCES workflow, demonstrating configuration file inputs and descriptions of processes and outputs for each PISCES subcommand (index, submit, summarize-expression and summarize-qc).

   digraph {
    rankdir=TB;
    abundance [label="<f0> TMM(TPM)|<f1> counts|<f2> transcript metrics", shape=record];
    de [label="<f0> log2(fold-change)|<f1> p-value", shape=record];
    qc [label="<f0> nucleotide metrics|<f1> genotype metrics", shape=record];
    index [label="<f0> salmon index|<f1> gene-transcript mapping|<f2> transcript/intron/intergene models", shape=record];
    "config.json" -> "pisces index" [arrowhead=none];
    "metadata.csv" -> "pisces submit" [arrowhead=none];
    "metadata.csv" -> "pisces summarize-expression" [arrowhead=none];
    "metadata.csv" -> "pisces summarize-qc" [arrowhead=none];
    "contrasts.csv" -> "pisces summarize-expression" [arrowhead=none];
    "pisces submit" -> "pisces summarize-expression";
    "pisces submit" -> "pisces summarize-qc";
    "pisces summarize-expression" -> abundance;
    "pisces summarize-expression" -> de;
    "pisces summarize-qc" -> qc;
    "pisces index" -> index;
    index -> "pisces submit";
    }

The PISCES workflow (:numref:`workflow_graph`) is designed to be **flexible**,
requiring minimal setup to describe an experiment; **scalable**, running on one to many samples either
locally or on a compute cluster; **fast**, using alignment-free methods
to generate QC metrics; and **reproducible**, driven by easily generated
configuration files and versioned transcriptome annotation, with automated
generation of transcriptome index files. We demonstrate each step of the PISCES workflow on a real experiment
(:numref:`treatment_graph`) and show examples of PISCES outputs on real data.

Results
=======

As a demonstration of the speed and simplicity of the PISCES workflow we
reprocess an NCBI SRA study consisting of two breast cancer cell lines, MCF7 and
T47D, with Y537S or D538G mutations introduced in the ESR1 (estrogen receptor)
gene and treated with Vehicle or Estrogen (:numref:`treatment_graph`).

.. _treatment_graph:
.. graphviz::
  :caption: Treatment information for SRP093386

   digraph {
    rankdir=LR;
    compound=true;
    subgraph clusterbiospecimen {
      label="CellLine";
      "MCF7";
      "T47D";
    }
    subgraph clustergenotype {
      label="Genotype";
      geno [label="<f0> Y537S|<f1> D538G|<f2> WT", shape=record];
      }
    subgraph clustertreatment {
      label="Treatment";
      "Vehicle" [shape=diamond];
      "Estrogen" [shape=diamond];
      }
    "MCF7" -> "geno" [ltail=clusterbiospecimen, lhead=clustergenotype, arrowhead=none];
    "geno" -> "Vehicle" [ltail=clustergenotype, lhead=clustertreatment, arrowhead=none];
    "T47D" -> "geno" [ltail=clusterbiospecimen, lhead=clustergenotype, arrowhead=none];
    "geno" -> "Estrogen" [ltail=clustergenotype, lhead=clustertreatment, arrowhead=none];
   }

``pisces index`` builds custom transcriptome indices
----------------------------------------------------

PISCES is distributed with configuration for human, mouse, and human/mouse
chimera transcriptomes (for xenograft experiments) derived from Gencode v32 and
Gencode vM23 :cite:`harrow2006gencode`. Custom transcriptomes can easily be defined in a JSON
configuration file format. Extra sequences (such as bacterial genomes, spike ins for library
prep, or genetic knock-ins in the experimental design) can be added to
transcriptomes as arbitrary FASTA formatted files. PISCES builds transcript
sequences using the GTF file format by first creating all defined transcript
models per gene, and then concatenating all unique intron sequences to form one
intronic transcript model per gene. In this way PISCES enables simultaneous
quantification of both fully processed mRNA and contaminating pre-mRNA
transcripts. Intergenic regions can also be captured, allowing quantification of
all transcriptomic compartments of an organism’s genomic assembly. PISCES
leverages the pufferfish index :cite:`Almodaresi191874` enhancements in salmon
version 1.0 and above which allows indexing of vastly larger transcript
sequences.

CSV formatted metadata files are input for PISCES
-------------------------------------------------

We include a metadata file describing libraries generated from MCF7 cell lines
in :cite:`bahreini2017mutation` as a test script, distributed with PISCES. The
following analysis, reproduced from this test script, serves as an example of a
typical RNAseq workflow from FASTQ files to differential expression calls. In
this example, PISCES uses the NCBI SRA toolkit to directly obtain sequencing
data from NCBI servers, although paths to local data can be substituted in
*Fastq1*, and optionally *Fastq2* columns.

.. csv-table::
  :name: sra_table
  :file: analysis/metadata_exp.csv
  :header-rows: 1

:numref:`sra_table`: Experimental metadata such as treatment variables and
covariates of interest (cell line, time point etc.) are described in CSV format.

.. csv-table::
  :name: contrasts_table
  :file: analysis/contrasts.csv

:numref:`contrasts_table`: Contrasts for differential expression analysis are
defined using the variables present in the metadata table. Columns in this file do not have headers, but 
correspond to ``covariate name``, ``experimental value``, ``control value``.

With the above two files, an example PISCES workflow is only a few commands, and
can be run easily on most Linux and MacOS machines.

.. code-block:: shell

  $ pip install novartis-pisces
  $ pisces index
  $ pisces submit -m metadata.csv
  $ pisces summarize-qc -m metadata.csv -f fingerprint.txt
  $ pisces summarize-expression -m metadata.csv -f contrasts.csv -d "~Treatment + Genotype"
  
  
A formula for linear modeling of the experiment (**-d**) specifies the
**Treatment** and **Genotype** covariates, and the interaction between the two,
and p-values and log2 fold change are computed from the fit for each comparison
defined in the contrasts CSV. The resulting summarized expression and QC data is
output to tab separated text files 
suitable for further analysis and visualization using a suitable analytics
environment such as Jupyter or RStudio.

``pisces submit`` runs PISCES on a DRMAA enabled HPC
----------------------------------------------------

The PISCES ``pisces submit`` command simplifies the task of executing the main PISCES
workflow on a typical compute cluster, meaning a grid compute cluster using SGE,
UGE, SLURM, Torque, or another DRMAA enabled job scheduler. The ``pisces submit``
command takes an appropriately formatted metadata file similar to :numref:`sra_table` as input and
dispatches multiple ``pisces run`` jobs, passing along any extra arguments a user
supplies, as well as any cluster resource limits such as maximum job runtime,
maximum memory, and number of CPU cores requested. `pisces submit` then monitors
job log output and reports job status (number of jobs queued, running,
finished), enabling a user to quickly process and track a large number of RNAseq
libraries in parallel.

``pisces summarize-qc`` aggregates genetic QC metrics
-----------------------------------------------------

Sample quality and identity is paramount for interpretable RNAseq results, but
is missing from most alignment-free workflows. Here, read sequence quality
metrics are assessed using fastqp version 0.3.4. For studies involving human
samples, PISCES uses exact kmer counting (k=21) to generate a VCF per sample
based on 226 chosen high minor-allele frequency SNPs. Log-odds scores
representing the likelihood that two samples are genetically identical are used
to generate a sample identity matrix. To demonstrate PISCES’ ability to match
genetic identity even across independent experiments, we included breast cancer
cell lines from the Cancer Cell Line Encyclopedia [12] with our example SRA
analysis. In (:numref:`fingerprint`) we show the genetic fingerprints of the
CCLE :cite:`ghandi2019next` T47D and MCF7 cell lines match the expected clusters
from the SRA experiment.  Additionally, one can detect cell line pairs that were
derived from the same patient (the KPL1 cell line which is a derivative of MCF7,
and the autologous pair AU565 and SKBR3.)

.. figure:: analysis/fingerprint_plot.png
  :name: fingerprint

  Sample identity matrix showing two clusters of samples, corresponding to
libraries derived from MCF7 and T47D cell line models. The p_same score is a
log-odds value where 0 indicates an exact genetic match and more negative values
indicate a lower probability of genetic identity. Breast cancer cell lines from
the publicly available CCLE are shown, demonstrating PISCES ability to integrate
datasets across experiments. 

``pisces summarize-expression`` aggregates transcript and gene abundance estimates
----------------------------------------------------------------------------------

PISCES uses salmon :cite:`Patro2017` for estimating transcript and gene
abundances from RNAseq libraries. 
Salmon is computationally efficient as well as accurate in assignment of reads
to transcripts :cite:`srivastava2020alignment`. PISCES includes Salmon version
1.3.0. Library
type parameters (strand-specific, read pairs, and read pair orientation) are
inferred from input FASTQ files. Salmon is run with default parameters as well
as specifying 
``--seqBias``, ``--gcBias``, and ``--useVBOpt``. The ``pisces
summarize-expression`` subcommand 
aggregates transcriptomic counts and TPM estimates for multiple samples and
summarizes these 
values to gene-level data matrices with appropriate column names defined in the
metadata CSV file. 
TMM normalization is applied to the TPM estimates to control for highly
expressed outlier transcripts. 
If unprocessed transcript models are specified during ``pisces index``, data
matrices for 
"intronic" and "intergene" transcriptomic compartments are output and a summary
of percent coding, intergenic and 
intronic sequence estimates can be used to determine relative sample QC
thresholds :numref:`qc_plots_1`.

.. figure:: analysis/qc_plot.png
  :name: qc_plots_1
  
  Percent coding, intronic, and intergenic read content is calculated salmon counts of processed transcript models, 
  transcript models of introns, and intergenic sequences.
  
If contrasts (:numref:`contrasts_table`) and a formula are specified, differential gene expression is calculated and output 
to a tidy data table. These in turn can be visualized using standard graphing techniques, for example as a volcano plot :numref:`volcano_plot`.

.. figure:: analysis/volcano_plot.png
  :name: volcano_plot
  
  Volcano plot of differentially expressed genes (DEGs) identified by DESeq2, with **p < 0.01**.  
  DEGs previously identified from :cite:`bahreini2017mutation` shown with the fold change stated in 
  the original publication. 
  
  
``pisces`` genomic indexing and repeat masking
----------------------------------------------

During development of PISCES it became clear that the expression of certain
genes was detected under biologically improbable contexts (such as aberrant immune
cell marker gene expression in non-immune cancer cell culture). Upon closer
inspection it appeared that many of these genes with unexpectedly high
expression were mapping reads to repetitive regions of the transcript sequence.
These repetitive regions are masked by RepeatMasker in the UCSC repeat masked genome assemblies,
and so by default the pisces index command will hard-mask any soft-masked
characters in the input genome assembly. The effect of repeat masking on salmon
TPM estimates is that many genes with apparently low-to-medium abundance are
significantly reduced, indicating their TPM estimates were driven primarily by
repeat sequences and were likely a result of incorrectly assigned reads 
(:numref:`repeatmask_plot`).

.. figure:: analysis/repeatmask_plot.png
  :name: repeatmask_plot
  
  Median estimated read counts per gene for MCF7 models (:numref:`sra_table`) for either 
  exon-centric transcript models (x-axis) or exon-centric transcript models with repeat-masking (y-axis) demonstrates 
  inflation (right shifted bins) of gene counts in the exon-centric data. Median counts were 
  binned in quarter-log increments with a pseudo-count added to every observation to avoid 
  undefined values.
  
Another way to visualize the effect of repeat masking is to examine its impact on differential 
gene expression results. In :numref:`de_plot_nomask` this effect is demonstrated between  exon-centric 
transcriptome indices with and without repeat masking. Off-axis genes contain reads that 
are potentially mis-mapped due to repetitive sequence content. Copies of these repetitive sequences 
may be found throughout the genome assembly, and without repeat masking, reads that map 
well to these sequences may be incorrectly assigned.

.. figure:: analysis/de_plot_nomask.png
  :name: de_plot_nomask
  
  Differentially expressed genes identified in D538G or Y537S vs. Wild-type MCF7 models (:numref:`sra_table`) 
  for either exon-centric transcript models (x-axis) or exon-centric transcript models with repeat-masking (y-axis). Off-axis 
  genes indicate potential mis-mapped reads due to repetitive sequence. DEGs with p < 0.01 are shown. Genes 
  with a greater than 2 fold difference between methods are labeled. Points are colored by the ratio of median gene counts 
  for exon-centric vs exon-centric repeat-masked methods.

Surprisingly, we can achieve almost exactly the same effect as repeat-masking by simply including
intronic and intergenic sequences to build a genome-centric transcript index for salmon. 
This may be due to missing transcripts in our transcriptome model that contain
repeats; by including these repeats in our intronic and intergenic sequences
they are no longer assigned to the transcript that did have these repeats. LINE and SINE repeats 
are commonly expressed in several cancer types :cite:`burns2017transposable`, and transcript 
models that contain these transposable elements, especially in 5' and 3' untranslated regions (UTRs) 
could be susceptible to such read mis-mapping. :numref:`de_plot` shows that changes to differentially expressed genes are
minimal between the repeat-masked exon-centric indexing method and the
genome-centric transcriptome indexing method. By compartmentalizing the entire
genome assembly we may reduce the number of spurious mappings arising from
repetitive genomic sequences.

.. figure:: analysis/de_plot.png
  :name: de_plot
  
  Differentially expressed genes identified in D538G or Y537S vs. Wild-type MCF7 models (:numref:`sra_table`) 
  for either exon-centric transcript models with repeat-masking (x-axis) or genome-centric transcript models (y-axis). Off-axis 
  genes indicate potential mis-mapped reads due to repetitive sequence. DEGs with p < 0.01 are shown. Genes 
  with a greater than 2 fold difference between exon-centric (repeat-masked) and genome-centric methods are labeled. Points are colored by the ratio of median gene counts 
  for exon-centric repeat-masked vs genome-centric methods.
  
Based on these results we recommend the use of a genome-centric transcriptome index in PISCES.

Discussion
==========

The motivation for development of PISCES is to enable the efficient, reproducible, and
automated analysis of the majority of RNAseq experiments. To this end, the main
design goals were ease of installation and use, built-in
cluster submission to enable rapid job processing, clarity and reproducibility
through use of configuration files, and a design that encourages
users to describe experiments using metadata files that can be leveraged at
multiple steps in the PISCES workflow. In this way users are encouraged toward
reproducible description of their RNAseq analysis. 

Traditional alignment-based RNAseq pipelines trade longer run times for more
detailed alignment of each sequence read. The gains in computing efficiency
acheived by pseudoalignment based transcript quantification come at the cost of
this detailed alignment information. This required adaptation of existing RNAseq
QC methods that rely on alignment, to keep the total run time of the PISCES
pipeline close to the run time of its' slowest individual component. 
Integration of external datasets (e.g. TCGA, GTEx, or SRA projects) can be
complicated by varying quality control from different sites, data transfer
concerns due to disparate (local and remote) data sources, and limitations due
to compute times necessary to run traditional tools. PISCES enables large scale
integration of datasets, processing of high throughput RNAseq generated from
large plate format screens, fast quantification against multiple custom 
transcriptomes, and fast, reliable reprocessing of libraries whenever new
transcriptome annotation or genome assemblies become available.

In addition to gains in computational efficiency, PISCES enables rapid analysis
of routine RNAseq experiments, through automated modeling of differential gene
expression and generation of clean tables for visualizing gene expression,
sample clustering, genetic identity, and QC metrics. We demonstrate a novel genome-centric 
transcriptome indexing strategy that allows integrated transcript QC measures of intronic 
and intergenic sequences while also avoiding read mismapping that may occur with exon-centric transcriptome 
quantification. PISCES is available on GitHub_ at 
https://github.com/Novartis/pisces, and on the python package index (PyPI_) at https://pypi.org/project/novartis-pisces.


.. bibliography:: references.bib
  :style: unsrt
