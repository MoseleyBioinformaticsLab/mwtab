Release History
===============

2.0.0
~~~~~
-Can now read duplicate keys in "Additional sample data" and reproduce it in write, will validate it as an error.
-Writing out now ensures correct key ordering for JSON.
-Validation now validates the main sections not just the internals of them.
-Validate now checks that metabolites in the Data section are in the Metabolites section and vice versa.
-Batch processing from the command line is more fault tolerant and won't stop the batch for 1 bad file.
-Improved tokenizer so more files can be read in without error.
-Changed schema validation to use jsonschema instead of schema.
-Added validations for METABOLITES columns that try to give warnings for bad values, for example 'kegg_id' column should all be something like C00000.
-Expanded the standard column name functionality to look for many more column names than in the previous version and do it in a much more robust way.
-Changed validation error messages to be aware of whether the input file is JSON or mwTab and print accordingly.
-Reduced many spurious validation messages.
-Added a validation for when certain columns are found in METABOLITES, to look for the implied pair and warn if it isn't there. For example, retention_index and retention_index_type.
-Added validations on some values, such as gender.
-Many more various minor validations were added.


1.2.5.post1 (2022-05-11)
~~~~~~~~~~~~~~~~~~~~~~~~

**Improvements**

- Add citation information to GitHub repository.

    - Adds CITATION.cff file with citation info.


1.2.5 (2022-03-18)
~~~~~~~~~~~~~~~~~~

**Improvements**

- Updates ``mwschema.py`` and ``validator.py`` modules to match Metabolomics Workbench's mwTab File Format
Specification Version 1.5 (March 2022).

    - Adds optional NMR_RESULTS_FILE field to NMR block.

    - Adds optional MS_COMMENTS field to MS block.

    - Removes requirement for there to be data results for every sample listed in the Study Design
(SUBJECT_SAMPLE_FACTORS). Allows for instances where samples have technical issues preventing data from being provided.


1.2.4 (2022-01-07)
~~~~~~~~~~~~~~~~~~

**Improvements**

- Adds check for blank source files when parsing to create ``mwtab`` objects.


1.2.3 (2021-11-02)
~~~~~~~~~~~~~~~~~~

**Bugfixes**

- Removes hard coding of version number in ``validator.validate_file()`` method.

- Removes mention of Python 3.4 support in README.


1.2.2 (2021-10-22)
~~~~~~~~~~~~~~~~~~

**Improvements**

- Migrates Continuous Integration (CI) from Travis CI to GitHub Actions.

    - Adds ``.github/workflows/`` folder which contains .yml files for workflows.

        - Adds ``build.yml`` to folder for testing build with pytest.

        - Adds ``codecov.yml`` to folder for generating/uploading code coverage info to codecov.io
          (https://app.codecov.io/gh/MoseleyBioinformaticsLab/mwtab).

    - Changes build and codecov badges to match new sources.


1.2.1 (2021-09-03)
~~~~~~~~~~~~~~~~~~

**Improvements**

- Updates format of ``~mwtab.mwtab.validate_file()`` validation log generated during validation.

    - Includes metadata header in validation logs containing; datetime, mwtab version, file source, study id, analysis
      id, and file format.

    - Minor changes to error messages for MS(NMR)_METABOLITE_DATA, NMR_BINNED_DATA, and SUBJECT_SAMPLE_FACTORS sections.

**Bugfixes**

- Fixes error where pytests for ``~mwtab.mwtab.validate_file()`` method were repeatedly using the same text files for
validation rather than both the test text and JSON files.

- Verbose file validation enabled in commandline.

- Default value given to ``base_url`` parameter in ``~mwtab.mwatb._pull_study_analysis()`` methods.


1.0.1 (2021-03-06)
~~~~~~~~~~~~~~~~~~

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
