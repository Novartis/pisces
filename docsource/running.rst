Running
=======
Novartis HPC environment
------------------------
PISCES is packaged as a ``module`` file which must be added to your ``$MODULEPATH``. You can set this in your ``${HOME}/.bashrc`` file for persistence:

.. code:: shell

  export MODULEPATH=${MODULEPATH}:/usr/prog/modules/all:/cm/shared/modulefiles/:/usr/prog/onc/modules/

You can then load PISCES by running:

.. code:: shell

  $ module load pisces/2018.04.1

  or

  $ module load pisces # for the latest version

If you'll be submitting jobs to the HPC compute notes you'll also need ``uge``:

.. code:: shell

  $ module load uge


PISCES index files are pre-built on the NIBR HPC environment, so you can directly skip to the section on `pisces submit`_.

.. _index_example:

pisces index
------------
.. code:: shell

    $ pisces index

Builds the default PISCES :doc:`config`.
This will build salmon index files for human, mouse, and human-mouse sample types.

.. note::

	If index output folders exist, they will not be overwritten, as PISCES assumes the index has already been built.

You can also pass in a custom ``--config`` file:

.. code:: shell

    $ pisces --config index

.. _run_example:

pisces run
----------
.. code:: shell

    $ pisces run -fq1 lane1_R1_001.fastq.gz lane1_R1_002.fastq.gz \
                 -fq2 lane1_R2_001.fastq.gz lane1_R2_002.fastq.gz

In the most basic form, you can specify only the fastq files (as a list
of forward and reverse reads) and other parameters will be auto-detected
or selected from default values. Either paired or unpaired libraries are
allowed. If the data are unpaired, just pass fastq files using ``-fq1``.

.. todo::

  Fix the link to default config file

Data and program paths are defined using a default 
:doc:`config` which can be specified at runtime using the
``--config`` argument.

.. code:: shell

    $ pisces --config config.json run \
             -fq1 lane1_R1_001.fastq.gz lane1_R1_002.fastq.gz \
             -fq2 lane1_R2_001.fastq.gz lane1_R2_002.fastq.gz

Sample name (``-n, --name``), output directory (``-o, --out``) and total
number of CPU threads to utilize (``-p, --threads``) may be specified
explicitly, or default to automatic values.

.. code:: shell

    $ pisces run -fq1 lane1_R1_001.fastq.gz lane1_R1_002.fastq.gz \
                 -fq2 lane1_R2_001.fastq.gz lane1_R2_002.fastq.gz \
                 -p 12 \
                 -o PISCES_output_sample1

.. command-output:: pisces run --help

.. _summarize_example:

pisces summarize-expression
---------------------------

.. code:: shell

    $ pisces summarize Sample1/PISCES Sample2/PISCES Sample3/PISCES ...

or

.. code:: shell

    $ pisces summarize -m metadata.csv

You can summarize transcript-level expression to gene-level and make TPM
and counts matrices using ``pisces summarize``. Required arguments are
the directories specified as ``--out`` from ``pisces run``. Optionally
you can supply a metadata matrix in CSV format similar to `this
example <data/metadata_example.csv>`__:

+------------+-------------+-------------+
| SampleID   | Treatment   | Timepoint   |
+============+=============+=============+
| Sample1    | DMSO        | 1h          |
+------------+-------------+-------------+
| Sample2    | DMSO        | 1h          |
+------------+-------------+-------------+
| Sample3    | DMSO        | 1h          |
+------------+-------------+-------------+
| Sample4    | Dox         | 1h          |
+------------+-------------+-------------+
| Sample5    | Dox         | 1h          |
+------------+-------------+-------------+
| Sample6    | Dox         | 1h          |
+------------+-------------+-------------+
| Sample7    | DMSO        | 4h          |
+------------+-------------+-------------+
| Sample8    | DMSO        | 4h          |
+------------+-------------+-------------+
| Sample9    | DMSO        | 4h          |
+------------+-------------+-------------+
| Sample10   | Dox         | 4h          |
+------------+-------------+-------------+
| Sample11   | Dox         | 4h          |
+------------+-------------+-------------+
| Sample12   | Dox         | 4h          |
+------------+-------------+-------------+

When supplying a ``--metadata`` file you can specify the ``--group-by``
option to group samples (e.g. Timepoint) before normalizing using the
``--norm-by`` variable (e.g. Treatment) with the ``--control-factor``
(e.g. DMSO) as the set of control samples to normalize to. You can also
pass a formula for differential expression using DESeq2 by specifying
``--deseq-formula`` such as
``--deseq-formula "~ Treatment + Treatment:Timepoint"``. The
``--spotfire-template`` option copies a template Spotfire file useful
for visualizing the resulting data matrices.

By default ``pisces summarize`` matches metadata to input sample
directories based on the order of directories passed as positional
arguments. E.g:
``pisces summarize -m metadata.csv /Sample1 /Sample2 ...``. Sometimes
this is cumbersome, so there are two options for encoding input
locations in the metadata file:

As paths to ``pisces run`` output directories:

+------------+-------------+-------------------------+
| SampleID   | Treatment   | Directory               |
+============+=============+=========================+
| Sample1    | DMSO        | /path/to/PISCES\_run1   |
+------------+-------------+-------------------------+
| Sample2    | DMSO        | /path/to/PISCES\_run2   |
+------------+-------------+-------------------------+

As paths to salmon "quant.sf" files:

+------------+-------------+----------------------------------+
| SampleID   | Treatment   | QuantFilePath                    |
+============+=============+==================================+
| Sample1    | DMSO        | /path/to/PISCES\_run1/quant.sf   |
+------------+-------------+----------------------------------+
| Sample2    | DMSO        | /path/to/PISCES\_run2/quant.sf   |
+------------+-------------+----------------------------------+

.. command-output:: pisces summarize-expression --help

.. _qc_example:

pisces summarize-qc
-------------------

QC tables are created using the ``pisces qc`` command. PISCES samples
are discovered recursively for each directory passed to the tool.

.. code:: shell

    $ pisces qc . \
                --spotfire-template QC.dxp \
                --tab QC.table.txt \
                --tall QC.skinny.txt \
                --fingerprint fingerprint_identities.txt

or

.. code:: shell

    $ pisces qc --metadata metadata.csv \
                --spotfire-template QC.dxp \
                --tab QC.table.txt \
                --tall QC.skinny.txt \
                --fingerprint fingerprint_identities.txt

Note that directories are searched recursively and so it is sufficient
to pass in the top level directory when all PISCES runs in the directory
are desired.

.. command-output:: pisces summarize-qc --help

.. _submit_example:

pisces submit
------------------
PISCES contains a command for running multiple ``pisces run`` jobs on a DRMAA-aware
compute cluster (sge, uge, slurm). Jobs are specified using the ``metadata.csv`` table
by adding data locations for the FASTQ files. Extra arguments to ``pisces run`` are passed to
``pisces submit`` and appended to each job before submission to the cluster. The DRMMA library
needs to be accessible in your environment: ``export DRMAA_LIBRARY_PATH=/cm/shared/apps/uge/current/lib/lx-amd64/libdrmaa.so``.

.. code:: shell

    $ pisces submit --metadata metadata.csv [pisces run args]

.. raw:: html

    <script type="text/javascript" src="https://asciinema.org/a/mdA7nJ91l8BjIeFjNmmWBjds4.js" id="asciicast-mdA7nJ91l8BjIeFjNmmWBjds4" async></script>

After job submission, ``pisces submit`` will monitor the progress of submitted
jobs. If you want to exit this command, pressing ``Ctrl+C`` will prompt whether
to delete the current jobs. Job progress (running, completion, or failure) can
be checked at any time by re-running ``pisces submit`` in the directory where
``pisces submit`` was originally run. If you need to later re-run ``pisces submit`` in
the same directory you must first delete the ``.pisces`` directory.

.. command-output:: pisces submit --help
