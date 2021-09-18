import pytest
import mwtab


@pytest.mark.parametrize("files_source", [
    "tests/example_data/mwtab_files/ST000122_AN000204.json",
    "tests/example_data/mwtab_files/ST000122_AN000204.txt"
])
def test_validate(files_source):
    """Test method for validating passing mwTab and JSON files from Metabolomics Workbench.
    :param files_source: File path to Metabolomics Workbench file to be validated.
    :type files_source: :py:class:`str` or
    """
    mwfile = next(mwtab.read_files(files_source))
    _, validation_log = mwtab.validate_file(mwfile, metabolites=False)
    assert len(validation_log.split('\n')) == 9


@pytest.mark.parametrize("file_source", [
    "tests/example_data/validation_files/ST000122_AN000204_error_1.txt",
    "tests/example_data/validation_files/ST000122_AN000204_error_1.json"
])
def test_validate_subject_sample_factors(file_source):
    mwfile = next(mwtab.read_files(file_source))
    _, validation_log = mwtab.validate_file(mwfile, metabolites=False)
    assert "missing Subject ID" in validation_log
    assert "missing Sample ID" in validation_log
    assert "missing value for Factor" in validation_log


@pytest.mark.parametrize("file_source", [
    "tests/example_data/validation_files/ST000122_AN000204_error_2.txt",
    "tests/example_data/validation_files/ST000122_AN000204_error_2.json"
])
def test_validate_subject_sample_factors(file_source):
    mwfile = next(mwtab.read_files(file_source))
    _, validation_log = mwtab.validate_file(mwfile, metabolites=False)
    assert "Section missing data entry for sample(s):" in validation_log
    assert "SUBJECT_SAMPLE_FACTORS: Section missing sample ID(s)" in validation_log


@pytest.mark.parametrize("file_source", [
    "tests/example_data/validation_files/ST000122_AN000204_error_3.txt",
    "tests/example_data/validation_files/ST000122_AN000204_error_3.json"
])
def test_validate_metabolites(file_source):
    mwfile = next(mwtab.read_files(file_source))
    _, validation_log = mwtab.validate_file(mwfile)
    assert "which matches a commonly used field name" in validation_log


@pytest.mark.parametrize("file_source", [
    "tests/example_data/validation_files/ST000122_AN000204_error_4.txt",
    "tests/example_data/validation_files/ST000122_AN000204_error_4.json"
])
def test_validate_schema(file_source):
    mwfile = next(mwtab.read_files(file_source))
    _, validation_log = mwtab.validate_file(mwfile)
    assert "does not match the allowed schema" in validation_log


@pytest.mark.parametrize("file_source", [
    "tests/example_data/mwtab_files/ST000122_AN000204.json"
])
def test_validation_log_local(file_source):
    mwfile = next(mwtab.read_files(file_source))
    _, validation_log = mwtab.validate_file(mwfile)
    # assert "mwtab version: {}".format(mwtab.__version__) in validation_log
    assert "Source:        {}".format(file_source) in validation_log
    assert "Study ID:      {}".format("ST000122") in validation_log
    assert "Analysis ID:   {}".format("AN000204") in validation_log
    assert "File format:   {}".format("json") in validation_log


@pytest.mark.parametrize("file_source", [
    "2"
])
def test_validation_log_web(file_source):
    mwfile = next(mwtab.read_files(file_source))
    _, validation_log = mwtab.validate_file(mwfile, metabolites=False)
    # assert "mwtab version: {}".format(mwtab.__version__) in validation_log
    assert "Source:        {}".format("https://www.metabolomicsworkbench.org/rest/study/analysis_id/AN000002/mwtab/txt")\
           in validation_log
    assert "Study ID:      {}".format("ST000002") in validation_log
    assert "Analysis ID:   {}".format("AN000002") in validation_log
    assert "File format:   {}".format("txt") in validation_log