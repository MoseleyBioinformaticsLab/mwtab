import pytest
import mwtab


@pytest.mark.parametrize("files_source", [
    "tests/example_data/mwtab_files/ST000122_AN000204.json",
    "tests/example_data/mwtab_files/ST000122_AN000204.txt"
])
def test_validate(files_source):
    mwfile = next(mwtab.read_files(files_source))
    _, validation_errors = mwtab.validate_file(mwfile, metabolites=False)

    assert not validation_errors


@pytest.mark.parametrize("file_source", [
    "tests/example_data/validation_files/ST000122_AN000204_error_1.txt",
    "tests/example_data/validation_files/ST000122_AN000204_error_1.txt"
])
def test_validate_subject_sample_factors(file_source):
    mwfile = next(mwtab.read_files(file_source))
    _, validation_errors = mwtab.validate_file(mwfile, metabolites=False)

    assert "missing Subject ID" in validation_errors
    assert "missing Sample ID" in validation_errors
    assert "missing value for Factor" in validation_errors


@pytest.mark.parametrize("file_source", [
    "tests/example_data/validation_files/ST000122_AN000204_error_2.txt",
    "tests/example_data/validation_files/ST000122_AN000204_error_2.txt"
])
def test_validate_data_section(file_source):
    mwfile = next(mwtab.read_files(file_source))
    _, validation_errors = mwtab.validate_file(mwfile, metabolites=False)

    assert "\" missing data entry for" in validation_errors


@pytest.mark.parametrize("file_source", [
    "tests/example_data/validation_files/ST000122_AN000204_error_3.txt",
    "tests/example_data/validation_files/ST000122_AN000204_error_3.txt"
])
def test_validate_metabolites(file_source):
    mwfile = next(mwtab.read_files(file_source))
    _, validation_errors = mwtab.validate_file(mwfile)

    assert "which matches a commonly used field name" in validation_errors


@pytest.mark.parametrize("file_source", [
    "tests/example_data/validation_files/ST000122_AN000204_error_4.txt",
    "tests/example_data/validation_files/ST000122_AN000204_error_4.txt"
])
def test_validate_schema(file_source):
    mwfile = next(mwtab.read_files(file_source))
    _, validation_errors = mwtab.validate_file(mwfile)

    assert "does not match the allowed schema" in validation_errors
