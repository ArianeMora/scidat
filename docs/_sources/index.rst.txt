*******
sci-dat
*******

Sci-Download-Annotate-TCGA is a wrapper around the functions provided by TCGA and the `GDC data portal <https://portal.gdc.cancer.gov/>`_.
Long story short, I was needing to merge many of the data (RNAseq and DNA methylation) together from TCGA and I wanted
to keep track of the demographics of the patients to ensure I had a balanced dataset. I also wanted to easily find
genes in groups of patients with mutations. I found no easy ways to do these things, so I made this wrapper to be able to:

1) Create a dataframe of many RNAseq datasets from TCGA (and automatically download these)
2) Merge RNAseq and DNA methylation datasets so for each gene I could see a cross mode profile
3) Annotate each experiment with demographic information
4) Anotate each gene with mutation information and search for genes with specific mutations through the API.

This package provides the above in `python notebooks`, `R markdown`, and a `CLI`.

It is available under the `GNU General Public License (Version 3) <https://www.gnu.org/licenses/gpl-3.0.en.html>`_.

Please post questions and issues related to sci-dat on the
`Issues <https://github.com/ArianeMora/scidat/issues>`_  section of the GitHub repository.


Running sci-dat
===============

1. Install sci-dat (:ref:`Installing <installing>`)

2. View examples in Python

3. View examples in R

4. Look at CLI examples

Extending sci-dat
=================

1. Make a pull request on github.


Citing sci-dat
==================
Sci-dat can be cited as in :ref:`references`, where we also provide citations for the used tools (e.g. numpy and GDC-data portal).

.. toctree::
   :caption: Getting started
   :maxdepth: 1

   about
   installing/index


.. toctree::
   :caption: Running sci-dat
   :maxdepth: 1

   examples/examples
   examples/cli
   examples/notebook


.. toctree::
   :caption: About
   :maxdepth: 1

   faq
   changelog
   references
