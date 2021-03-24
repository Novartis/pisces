[![CI](https://github.com/Novartis/pisces/workflows/CI/badge.svg)](https://github.com/Novartis/pisces/actions?query=workflow%3ACI)
[![Sphinx Docs](https://github.com/Novartis/pisces/workflows/Sphinx%20Docs/badge.svg)](https://opensource.nibr.com/pisces/)
[![Package Builds](https://github.com/Novartis/pisces/workflows/Package%20Builds/badge.svg)](https://github.com/Novartis/pisces/actions?query=workflow%3A%22Package+Builds%22)
[![PyPI](https://badge.fury.io/py/novartis-pisces.svg)](https://pypi.python.org/pypi/novartis-pisces)

# PISCES: a package for quantitation and QC of large scale mRNA-seq datasets

PISCES is a pipeline for rapid transcript quantitation, genetic fingerprinting, and quality control assessment of RNAseq libraries using Salmon.

See the current documentation at https://opensource.nibr.com/pisces/. Read the manuscript at https://www.biorxiv.org/content/10.1101/2020.12.01.390575v1.

Matthew D Shirley, Viveksagar K Radhakrishna, Javad Golji, Joshua M Korn. (Preprint) PISCES: a package for rapid quantitation and quality control of large scale mRNA-seq datasets. bioRxiv. 2020.

## Quickstart

Installation:
```
$ pip install --user novartis-pisces
$ pisces_dependencies
```

Submitting jobs to an HPC cluster:
```
$ pisces submit -m metadata.csv
```

Summarizing QC metrics and gene expression:
```
$ pisces summarize-qc -m metadata.csv -f fingerprint.txt
$ pisces summarize-expression -m metadata.csv
# optionally specify contrasts and formula for DESeq2 differential gene expression
$ pisces summarize-expression -m metadata.csv -f contrasts.csv -d "~covariate1 + covartiate2"
```

Please submit bugs or questions to our issue tracker.
