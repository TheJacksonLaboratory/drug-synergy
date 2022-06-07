Welcome to DECNEO documentation!
================================

.. image:: https://badgen.net/badge/Open%20Source/Yes/blue?icon=github
    :target: https://github.com/sdomanskyi/decneo/
.. image:: https://img.shields.io/github/release/sdomanskyi/decneo.svg
    :target: https://github.com/sdomanskyi/decneo/releases
.. image:: https://badge.fury.io/py/decneo.svg
    :target: https://pypi.org/project/decneo
.. image:: https://readthedocs.org/projects/decneo/badge/?version=latest
    :target: https://decneo.readthedocs.io/en/latest/?badge=latest
.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.4417915.svg
   :target: https://doi.org/10.5281/zenodo.4417915
.. image:: https://img.shields.io/twitter/url?style=social&url=https%3A%2F%2Fdoi.org%2F10.1101%2F2021.01.04.425317
    :target: https://twitter.com/intent/tweet?original_referer=https%3A%2F%2Fwww.biorxiv.org%2F&ref_src=twsrc%5Etfw&text=Comberons%20from%20single%20cell%20transcriptomics%20in%20endothelial%20cells&tw_p=tweetbutton&url=https%3A%2F%2Fwww.biorxiv.org%2Fcontent%2F10.1101%2F2021.01.04.425317v1

DECNEO (DEndrogram of Co-expression with Evolutionary Conservation and Network Enrichment Ordering) is a new software tool that identifies general and tissue-specific receptor comberons in endothelial cells and provides visualization and statistical analysis.

.. Note:: 

    DECNEO is introduced in:

    **Comberons from single cell transcriptomics in endothelial cells.** Sergii Domanskyi, Alex Hakansson, Michelle Meng, Joshua S Graff Zivin, Carlo Piermarocchi, Giovanni Paternostro,  Napoleone Ferrara, *in review* (2021). BioRxiv preprint: https://doi.org/10.1101/2021.01.04.425317

The pipeline shown in :ref:`overview` is implemented in our new software DECNEO. We have developed a main and an alternative implementation of DECNEO that are open-source and maintained at https://github.com/sdomanskyi/decneo. Releases are published on the Zenodo archive at https://doi.org/10.5281/zenodo.4417915. 

The main implementation can be installed from PyPI using the command ``pip install decneo``, see :ref:`getting-started`. DECNEO modules are relying on Python numpy, scipy, pandas, scikit-learn, matplotlib and other packages (see :ref:`overview` for dependencies). The main implementation can leverage on a multicore large-memory High-Performance Computing Systems for fast processing of large datasets. The software was optimized for use with the HPCC supercomputer at Michigan State University but also tested on different operating systems and hardware.

.. seealso:: The source code of the alternative implementation is deposited at https://github.com/sdomanskyi/decneo/validation.

Contents
========

.. toctree::
   :maxdepth: 4

   overview
   gettingStarted
   mainClass
   data
   examples
   outputData


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
