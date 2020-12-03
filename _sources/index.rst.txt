PISCES Documentation
====================

PISCES is a pipeline for rapid transcript quantitation, genetic fingerprinting,
and quality control assessment of RNAseq libraries
using `Salmon <http://github.com/COMBINE-lab/salmon>`__.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   install.rst
   data_formats.rst
   managing_metadata.rst
   running.rst
   examples.rst
   manuscript.rst


Quickstart
==========

- Install the python module (see :doc:`install`)
- Create transcriptome and genome index files (see :ref:`index_example`)
- Create a CSV file to define your experiment (see :doc:`managing_metadata`)
- Submit your jobs to a DRMAA aware high performance computing cluster (see :ref:`submit_example`)
- Summarize gene/transcript expression and QC metrics in analysis-ready tables using :ref:`summarize_example` and :ref:`qc_example`.

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

Publication
===========

The manuscript is available on `bioRxiv <https://www.biorxiv.org/content/10.1101/2020.12.01.390575v1>`__ and the latest version can be viewed in this documentation: :doc:`manuscript`.

Matthew D Shirley, Viveksagar K Radhakrishna, Javad Golji, Joshua M Korn. `PISCES: a package for rapid quantitation and quality control of large scale mRNA-seq datasets <https://www.biorxiv.org/content/10.1101/2020.12.01.390575v1>`__. bioRxiv. 2020.

Development
===========

Pull requests and issues are welcomed on the pisces `Github repository <http://github.com/Novartis/pisces>`__.
