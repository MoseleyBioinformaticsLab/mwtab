mwtab
=====

.. image:: https://img.shields.io/pypi/l/mwtab.svg
   :target: https://choosealicense.com/licenses/bsd-3-clause-clear/
   :alt: License information

.. image:: https://img.shields.io/pypi/v/mwtab.svg
   :target: https://pypi.org/project/mwtab
   :alt: Current library version

.. image:: https://img.shields.io/pypi/pyversions/mwtab.svg
   :target: https://pypi.org/project/mwtab
   :alt: Supported Python versions

.. image:: https://readthedocs.org/projects/nmrstarlib/badge/?version=latest
   :target: http://mwtab.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation status

.. image:: https://github.com/MoseleyBioinformaticsLab/mwtab/actions/workflows/build.yml/badge.svg
   :target: https://github.com/MoseleyBioinformaticsLab/mwtab/actions/workflows/build.yml
   :alt: Build status

.. image:: https://codecov.io/gh/MoseleyBioinformaticsLab/mwtab/branch/master/graph/badge.svg?token=jhjMsP1qma
   :target: https://codecov.io/gh/MoseleyBioinformaticsLab/mwtab
   :alt: CodeCov

.. image:: https://img.shields.io/badge/DOI-10.3390%2Fmetabo11030163-blue.svg
   :target: https://doi.org/10.3390/metabo11030163
   :alt: Citation link

.. image:: https://img.shields.io/github/stars/MoseleyBioinformaticsLab/mwtab.svg?style=social&label=Star
   :target: https://github.com/MoseleyBioinformaticsLab/mwtab
   :alt: GitHub project

|

.. image:: https://raw.githubusercontent.com/MoseleyBioinformaticsLab/mwtab/master/docs/_static/images/mwtab_logo.png
   :width: 50%
   :align: center
   :target: http://mwtab.readthedocs.io/


The ``mwtab`` package is a Python library that facilitates reading and writing
files in ``mwTab`` format used by the `Metabolomics Workbench`_ for archival of
Mass Spectrometry (MS) and Nuclear Magnetic Resonance (NMR) experimental data.

The ``mwtab`` package provides facilities to convert ``mwTab`` formatted files into
their equivalent ``JSON`` ized representation and vice versa.  ``JSON`` stands for JavaScript
Object Notation, an open-standard format that uses human-readable text to transmit
data objects consisting of attribute-value pairs.

The ``mwtab`` package can be used in several ways:

   * As a library for accessing and manipulating data stored in ``mwTab`` format files.
   * As a command-line tool to convert between ``mwTab`` format and its equivalent
     ``JSON`` representation.


Citation
~~~~~~~~

When using ``mwtab`` package in published work, please cite the following papers:

   * Powell, Christian D., and Hunter NB Moseley. "The mwtab Python Library for RESTful
     Access and Enhanced Quality Control, Deposition, and Curation of the Metabolomics
     Workbench Data Repository." *Metabolites* 11.3 (2021): 163. doi:
     `10.3390/metabo11030163`_.

   * Smelter, Andrey and Hunter NB Moseley. "A Python library for FAIRer access and
     deposition to the Metabolomics Workbench Data Repository."
     *Metabolomics* 2018, 14(5): 64. doi: `10.1007/s11306-018-1356-6`_.


Links
~~~~~

   * mwtab @ GitHub_
   * mwtab @ PyPI_
   * Documentation @ ReadTheDocs_


Installation
~~~~~~~~~~~~

The ``mwtab`` package runs under Python 3.4+. Use pip_ to install.
Starting with Python 3.4, pip_ is included by default.


Install on Linux, Mac OS X
--------------------------

.. code:: bash

   python3 -m pip install mwtab


Install on Windows
------------------

.. code:: bash

   py -3 -m pip install mwtab


Upgrade on Linux, Mac OS X
--------------------------

.. code:: bash

   python3 -m pip install mwtab --upgrade


Upgrade on Windows
------------------

.. code:: bash

   py -3 -m pip install mwtab --upgrade


Quickstart
~~~~~~~~~~

.. code:: python

   >>> import mwtab
   >>>
   >>> # Here we use ANALYSIS_ID of file to fetch data from URL
   >>> for mwfile in mwtab.read_files("1", "2"):
   ...      print("STUDY_ID:", mwfile.study_id)
   ...      print("ANALYSIS_ID:", mwfile.analysis_id)
   ...      print("SOURCE:", mwfile.source)
   ...      print("Blocks:", list(mwfile.keys()))
   >>>


.. image:: https://raw.githubusercontent.com/MoseleyBioinformaticsLab/mwtab/master/docs/_static/images/mwtab_demo.gif
   :align: center


.. note:: Read the User Guide and the ``mwtab`` Tutorial on ReadTheDocs_
          to learn more and to see code examples on using the ``mwtab`` as a
          library and as a command-line tool.


License
~~~~~~~

This package is distributed under the BSD_ `license`.


.. _Metabolomics Workbench: http://www.metabolomicsworkbench.org
.. _GitHub: https://github.com/MoseleyBioinformaticsLab/mwtab
.. _ReadTheDocs: http://mwtab.readthedocs.io
.. _PyPI: https://pypi.org/project/mwtab
.. _pip: https://pip.pypa.io
.. _BSD: https://choosealicense.com/licenses/bsd-3-clause-clear/
.. _10.3390/metabo11030163: https://doi.org/10.3390/metabo11030163
.. _10.1007/s11306-018-1356-6: http://dx.doi.org/10.1007/s11306-018-1356-6
