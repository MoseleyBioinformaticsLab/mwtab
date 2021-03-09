.. :changelog:

Release History
===============


1.1.0 (2021-03-08)
~~~~~~~~~~~~~~~~~~~

**Bugfixes**

- Removes Python 2.7 testing/support and adds Python 3.9 testing/support.

    - Removes Python 2.7 Travis CI job in .travis.yml config file.
    - Removes Python 2.7 support listing in setup.py (for PyPI).
    - Adds Python 3.9 Travis CI job in .travis.yml config file.
    - Adds Python 3.9 support listing in setup.py (for PyPI).


1.0.1 (2021-03-06)
~~~~~~~~~~~~~~~~~~~

**Improvements**

- Updated ``~mwtab.mwtab.MWTabFile`` to match Metabolomics Workbench JSON
  format.

    - Internal dictionary representation now matches Metabolomics Workbench
      JSON format.
    - ``~mwtab.mwtab.MWTabFile.write()`` and
      ``~mwtab.mwtab.MWTabFile.write_str()`` methods now produce files
      consistent with Metabolomics Workbench's JSON format.

- Updated ``mwschema.py`` to be consistent with Metabolomics Workbench's
  updated `mwTab` format specification.

- Added ``mwrest.py`` module for working with Metabolomics Workbench's REST API.

    - Allows for additional data file to be requested through Metabolomics
      Workbench's REST API.

- Added ``mwextract.py`` module for extracting metadata and metabolites from
  `mwTab` formatted files.

- Updated ``validator.py``.

    - Validator now collects all present errors.
    - Performs detection of common field names in `#METABOLITES` blocks.

- Updated ``docs/tutorial.ipynb`` to document improved and updated package
  functionality.

- Updated `mwtab` package to include Python 3.8 support.


0.1.10 (2019-02-18)
~~~~~~~~~~~~~~~~~~~

**Bugfixes**

- Metabolomics Workbench started using HTTPS,
  update reading from ANALYSIS_ID to address the change.


0.1.9 (2018-04-21)
~~~~~~~~~~~~~~~~~~

**Improvements**

- Added citation link to `mwtab` package.


0.1.8 (2018-04-05)
~~~~~~~~~~~~~~~~~~

**Improvements**

- Added `mwtab` package logo.
- Minor update: Simplified section validation function.


0.1.7 (2017-12-07)
~~~~~~~~~~~~~~~~~~

**Improvements**

- Minor update: Included test for additional header line within `mwTab` files
  that may or may not be present.


0.1.4, 0.1.5, 0.1.6 (2017-11-13)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Improvements**

- Minor update: package README file examples. 
- Minor update: update README to properly render on PyPI.


0.1.3 (2017-09-14)
~~~~~~~~~~~~~~~~~~

**Bugfixes**

- Fixed bug in the command-line interface.
- Fixed bug in ``mwschema.py`` module definition causing validation to fail.
- Fixed validation optional argument (to ``read_files()`` generator) in order
  to validate mwTab formatted files before returning them.
- Fixed Python2/3 compatibility bug that uses ``bz2`` Python module.
- Fixed Python2/3 unicode/str compatibility bug in ``mwschema.py`` module.

**Improvements**

- Added Travis CI tests: https://travis-ci.org/MoseleyBioinformaticsLab/mwtab
- Added code coverage reports: https://codecov.io/gh/MoseleyBioinformaticsLab/mwtab


0.1.2 (2017-09-14)
~~~~~~~~~~~~~~~~~~

**Bugfixes**

- Fixed issue with mwTab formatted file printable representation.


0.1.1 (2017-09-12)
~~~~~~~~~~~~~~~~~~

**Improvements**

- Improved README display on PyPI.


0.1.0 (2017-09-12)
~~~~~~~~~~~~~~~~~~

- Initial public release.
