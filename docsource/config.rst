configuration file format
=========================

PISCES is driven by JSON configuration files. The key:value
pairs in this configuration file then are used by all subsequent ``pisces``
subcommands. Paths can be 
either local files, or HTTP/FTP URLs. Notice that chimeric organisms such as a human/mouse xenograft organisms are
easily defined by specifying each transcriptome as a separate GTF/FASTA pair. 
Example transcriptomes for human, mouse, and human/mouse xenograft are included
in the default PISCES distribution. If you need to support another organism, this
can easily be accomplished by adding your transcriptome of interest in a new
configuration file. The configuration file format follows:

.. note::

   PISCES will locate the `salmon` executable using the following precedence:

   1. `PISCES_SALMON` environment variable (must reference an executable binary)
   2. system `salmon` in `PATH`
   3. bundled `redist/salmon/bin/salmon`

   If none are found, execution will fail with `FileNotFoundError`.

.. literalinclude:: ../pisces/config.json