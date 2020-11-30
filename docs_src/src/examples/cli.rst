.. cli:

CLI
===

CLI instructions for running scidat.

The CLI allows users to directly to download and annotate information from TCGA. RNAseq and DNA methylation datasets
can be generated and joined using CLI. Genes with mutations can also be recorded.

Example:
--------
Here we show an example where we download RNAseq and DNA methylation data from a manifest file and then merging the
two datasets.


.. code-block:: bash

    scidat --md "data/" --mf "data/manifest.tsv" --cf "data/clinical.txt" --sf "data/sample_sheet.txt" --af "data/tcga_hsapiens_gene_ensembl-GRCh38.p13.csv" --gdc "data/./gdc-client" --dd "./" --o "./" --rna t --meth t --download t


Arguments
---------

.. argparse::
   :module: scidat
   :func: gen_parser
   :prog: scidat
   :nodefaultconst:
