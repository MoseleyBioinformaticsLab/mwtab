mwtab
=====

.. image:: https://img.shields.io/pypi/l/mwtab.svg
   :target: https://choosealicense.com/licenses/bsd-3-clause
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

.. image:: https://api.travis-ci.org/MoseleyBioinformaticsLab/mwtab.svg?branch=master
   :target: https://travis-ci.org/MoseleyBioinformaticsLab/mwtab
   :alt: Travis CI status

.. image:: https://codecov.io/gh/MoseleyBioinformaticsLab/mwtab/branch/master/graphs/badge.svg?branch=master
   :target: https://codecov.io/gh/MoseleyBioinformaticsLab/mwtab
   :alt: Code coverage information


The `mwtab` package is a Python library that facilitates reading and writing
files in `mwTab` format used by the `Metabolomics Workbench`_ for archival of
Mass Spectrometry (MS) and Nuclear Magnetic Resonance (NMR) experimental data.

The `mwtab` package provides facilities to convert `mwTab` formatted files into
their equivalent JSONized (JavaScript Object Notation, an open-standard format that
uses human-readable text to transmit data objects consisting of attribute-value pairs)
representation and vice versa.

The `mwtab` package can be used in several ways:

   * As a library for accessing and manipulating data stored in `mwTab` format files.
   * As a command-line tool to convert between `mwTab` format and its equivalent
     `JSON` representation.


Links
~~~~~

   * mwtab @ GitHub_
   * mwtab @ PyPI_
   * Documentation @ ReadTheDocs_


Installation
~~~~~~~~~~~~

The `mwtab` package runs under Python 2.7 and Python 3.4+,
use pip_ to install. Starting with Python 3.4, pip_ is included
by default.


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


.. note:: Read the User Guide and the `mwtab` Tutorial on ReadTheDocs_
          to learn more and to see code examples on using the `mwtab` as a
          library and as a command-line tool.


License
~~~~~~~

This package is distributed under the BSD_ `license`.


.. _GitHub: https://github.com/MoseleyBioinformaticsLab/mwtab
.. _ReadTheDocs: http://mwtab.readthedocs.io
.. _PyPI: https://pypi.org/project/mwtab
.. _pip: https://pip.pypa.io
.. _Metabolomics Workbench: http://www.metabolomicsworkbench.org
.. _BSD: https://choosealicense.com/licenses/bsd-3-clause
