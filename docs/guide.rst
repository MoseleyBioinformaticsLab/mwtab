User Guide
==========

Description
~~~~~~~~~~~

The ``mwtab`` package is a Python library that facilitates reading and writing
files in ``mwTab`` format used by the `Metabolomics Workbench`_ for archival of
Mass Spectrometry (MS) and Nuclear Magnetic Resonance (NMR) experimental data.

The ``mwtab`` package provides facilities to convert ``mwTab`` formatted files into
their equivalent JSONized (JavaScript Object Notation, an open-standard format that
uses human-readable text to transmit data objects consisting of attribute-value pairs)
representation and vice versa.

The ``mwtab`` package can be used in several ways:

   * As a library for accessing and manipulating data stored in ``mwTab`` format files.
   * As a command-line tool to convert between ``mwTab`` format and its equivalent
     ``JSON`` representation.

Installation
~~~~~~~~~~~~

The :mod:`mwtab` package runs under Python 2.7 and Python 3.4+.
Starting with Python 3.4, pip_ is included by default. To install
system-wide with pip_ run the following:

Install on Linux, Mac OS X
--------------------------

.. code:: bash

   python3 -m pip install mwtab

Install on Windows
------------------

.. code:: bash

   py -3 -m pip install mwtab

Install inside virtualenv
-------------------------

For an isolated install, you can run the same inside a virtualenv_.

.. code:: bash

   $ virtualenv -p /usr/bin/python3 venv  # create virtual environment, use python3 interpreter

   $ source venv/bin/activate             # activate virtual environment

   $ python3 -m pip install mwtab         # install mwtab as usual

   $ deactivate                           # if you are done working in the virtual environment

Get the source code
~~~~~~~~~~~~~~~~~~~

Code is available on GitHub: https://github.com/MoseleyBioinformaticsLab/mwtab

You can either clone the public repository:

.. code:: bash

   $ https://github.com/MoseleyBioinformaticsLab/mwtab.git

Or, download the tarball and/or zipball:

.. code:: bash

   $ curl -OL https://github.com/MoseleyBioinformaticsLab/mwtab/tarball/master

   $ curl -OL https://github.com/MoseleyBioinformaticsLab/mwtab/zipball/master

Once you have a copy of the source, you can embed it in your own Python package,
or install it into your system site-packages easily:

.. code:: bash

   $ python3 setup.py install

Dependencies
~~~~~~~~~~~~

The :mod:`mwtab` package depends on several Python libraries. The ``pip`` command
will install all dependencies automatically, but if you wish to install them manually,
run the following commands:

   * docopt_ for creating :mod:`mwtab` command-line interface.
      * To install docopt_ run the following:

        .. code:: bash

           python3 -m pip install docopt  # On Linux, Mac OS X
           py -3 -m pip install docopt    # On Windows

   * schema_ for validating functionality of ``mwTab`` files based on ``JSON`` schema.
      * To install the schema_ Python library run the following:

        .. code:: bash

           python3 -m pip install schema  # On Linux, Mac OS X
           py -3 -m pip install schema    # On Windows


Basic usage
~~~~~~~~~~~

The :mod:`mwtab` package can be used in several ways:

   * As a library for accessing and manipulating data stored in ``mwTab`` formatted files.

      * Create the :class:`~mwtab.mwtab.MWTabFile` generator function that will generate
        (yield) a single :class:`~mwtab.mwtab.MWTabFile` instance at a time.

      * Process each :class:`~mwtab.mwtab.MWTabFile` instance:

         * Process ``mwTab`` files in a for-loop, one file at a time.
         * Process as an iterator calling the :py:func:`next` built-in function.
         * Convert the generator into a :py:class:`list` of :class:`~mwtab.mwtab.MWTabFile` objects.

   * As a command-line tool:

      * Convert from ``mwTab`` file format into its equivalent ``JSON`` file format and vice versa.
      * Validate data stored in ``mwTab`` file based on schema definition.

.. note:: Read :doc:`tutorial` to learn more and see code examples on using the :mod:`mwtab`
          as a library and as a command-line tool.


.. _pip: https://pip.pypa.io/
.. _virtualenv: https://virtualenv.pypa.io/
.. _docopt: https://pypi.org/project/docopt/
.. _schema: https://pypi.org/project/schema/
.. _Metabolomics Workbench: http://www.metabolomicsworkbench.org/
