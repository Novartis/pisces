PISCES metadata file format
---------------------------

.. tip::
   Column order does not matter to PISCES, but column names must match the required 
   fields exactly.

.. _submit_metadata_example:

.. table:: "pisces submit" metadata file example

+------------+-------------------------+---------------------+--------------------+
| SampleID   | Directory               | Fastq1              | Fastq2             |
+============+=========================+=====================+====================+
| Sample1    | /path/to/output_dir     | s1_R1_001.fastq.gz  | s1_R2_001.fastq.gz |
+------------+-------------------------+---------------------+--------------------+
| Sample2    | /path/to/output_dir     | s2_R1_001.fastq.gz  | s2_R2_001.fastq.gz |
+------------+-------------------------+---------------------+--------------------+

``SampleID``` is the unique identifier used to construct output folders, and as an 
identifier in :ref:`summarize_example` data table column headers. The ``Directory`` path 
points to the top level directory where PISCES outputs for a sample should be created. This 
may be a relative or absolute path. ``Fastq1`` is a required field for fragment (single end) 
sequencing libraries, and ``Fastq2`` is required to analyze paired end libraries. ``Fastq1`` 
and ``Fastq2`` paths can specify multiple files using a semicolon (``;``) separator:

+------------------------------------------+---------------------------------------+
| Fastq1                                   | Fastq2                                |
+==========================================+=======================================+
| s1_R1_001.fastq.gz;s1_R1_002.fastq.gz    | s1_R2_001.fastq.gz;s1_R2_002.fastq.gz |
+------------------------------------------+---------------------------------------+
| s2_R1_001.fastq.gz;s2_R1_002.fastq.gz    | s2_R2_001.fastq.gz;s2_R2_002.fastq.gz |
+------------------------------------------+---------------------------------------+

If NCBI Sequence Read Archive (SRA) accessions are specified, these must be added as ``SRR`` "run" accessions 
in the ``SRA`` column. If only SRA experiments are specified, the ``Fastq1`` column is 
optional.

+------------+-------------------------+---------------------+--------------------+-------------+
| SampleID   | Directory               | Fastq1              | Fastq2             | SRA         |
+============+=========================+=====================+====================+=============+
| Sample1    | /path/to/output_dir     | s1_R1_001.fastq.gz  | s1_R2_001.fastq.gz |             |
+------------+-------------------------+---------------------+--------------------+-------------+
| Sample2    | /path/to/output_dir     | s2_R1_001.fastq.gz  | s2_R2_001.fastq.gz |             |
+------------+-------------------------+---------------------+--------------------+-------------+
| Sample3    | /path/to/output_dir     |                     |                    | SRR000001   |
+------------+-------------------------+---------------------+--------------------+-------------+
| Sample4    | /path/to/output_dir     |                     |                    | SRR000002   |
+------------+-------------------------+---------------------+--------------------+-------------+


A metadata table such as this can be constructed using a bash script:

.. code:: shell

    $ ls Sample*
    Sample1:
    s1_R1_001.fastq.gz
    s1_R2_001.fastq.gz

    Sample2:
    s2_R1_001.fastq.gz
    s2_R2_001.fastq.gz

.. code:: shell

    echo "SampleID,Directory,Fastq1,Fastq2" > metadata.csv
    for dir in Sample*
      do
        fq1=$(ls $dir/*_R1_* | tr '\n' ';' | sed 's/;$//')
        fq2=$(ls $dir/*_R2_* | tr '\n' ';' | sed 's/;$//')
        printf "$dir,$dir/PISCES,$fq1,$fq2\n"
      done >> metadata.csv

.. tip::
   You may also find it easy to construct the metadata table using a spreadsheet editor.

Including analysis variables as metadata
========================================

PISCES utilizes variables defined in the metadata file when using :ref:`summarize_example` to 
run differential expression analysis, and for producing normalized fold-change tables. Any 
columns added to the file can be used in downstream analysis.

+------------+----------------+--------------+-------------------------+---------------------+--------------------+
| SampleID   | Treatment      | Timepoint    | Directory               | Fastq1              | Fastq2             |
+============+================+==============+=========================+=====================+====================+
| Sample1    | DMSO           | 4hours       | /path/to/output_dir     | s1_R1_001.fastq.gz  | s1_R2_001.fastq.gz |
+------------+----------------+--------------+-------------------------+---------------------+--------------------+
| Sample2    | Estrogen       | 12hours      | /path/to/output_dir     | s2_R1_001.fastq.gz  | s2_R2_001.fastq.gz |
+------------+----------------+--------------+-------------------------+---------------------+--------------------+

.. tip::
   For differential expression analysis in :ref:`summarize_example` it's often handy to 
   create "replicate group" variables composed of one or more treatment variables, e.g: 
   ``Treatment_Timepoint``. 

Specifying NCBI SRA projects
============================

You can easily create a metadata file for PISCES from the NCBI SRA "runinfo" format.
For example:

.. code:: shell

  $ wget -O - 'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=SRP093386' | \
    sed -e '1 s_Run,_SRA,_' -e '1 s_SampleName,_SampleID,_'  -e '1 s_Sample,_Directory,_' > metadata.csv

.. program-output:: wget --quiet -O - 'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=SRP093386' | sed -e '1 s_Run,_SRA,_' -e '1 s_SampleName,_SampleID,_'  -e '1 s_Sample,_Directory,_'
  :shell:
  :ellipsis: 8
