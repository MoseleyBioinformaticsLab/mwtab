.. :changelog:

Release History
===============

0.1.11 (2020-01-xx)
~~~~~~~~~~~~~~~~~~~

**Improvements**

- Added methods for working with Metabolomics Workbench's REST API.


**Bugfixes**

- Updated ``mwschema.py`` to be consistent with Metabolomics Workbench's
  updated `mwTab` format specification.

     - Adds `PROJECT_ID` as optional schema field in `#METABOLOMICS WORKBENCH` block.
- Updated ``validator.py`` to be consistent Metabolomics Workbench's
  updated `mwTab` format specification.

     - Now accepts `samples` and `factors` in a files `#MS_METABOLITE_DATA` and
       `#NMR_BINNED_DATA` blocks to be subsets of `sample` and `factors`
       listed in the file's `#SAMPLE_FACTORS` block.
- Fixed bug preventing files with multiple additional key value pairs in
  `_RESULTS_FILE` lines.

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
