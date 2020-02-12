import pytest
import mwtab
from collections import OrderedDict


@pytest.mark.parametrize("files_source", [
    "tests/example_data/validation_files/AN000001.txt",
    "tests/example_data/validation_files/AN000041.txt"
])
def test_validate(files_source):
    mwfile = next(mwtab.read_files(files_source))
    validation_errors = mwtab.validate_file(mwfile)
    assert type(validation_errors) == OrderedDict


@pytest.mark.parametrize("file_source", [[
    "tests/example_data/validation_files/AN000001_error_1.txt",
    "tests/example_data/validation_files/AN000001_error_2.txt"
]])
def test_validate_ms_samples(file_source):
    print(file_source)
    mwfile = next(mwtab.read_files(file_source[0]))
    validation_errors = mwtab.validate_file(mwfile, validate_factors=False, validate_features=False,
                                            validate_schema=False, validate_data=False, test=True)
    assert len(validation_errors) == 1
    test_error = KeyError("Missing key `Samples` in `MS_METABOLITE_DATA` block.")
    assert type(validation_errors[0]) == type(test_error) and validation_errors[0].args == test_error.args

    mwfile = next(mwtab.read_files(file_source[1]))
    validation_errors = mwtab.validate_file(mwfile, validate_factors=False, validate_features=False,
                                            validate_schema=False, test=True)
    assert len(validation_errors) == 3
    test_error = ValueError("Sample with no Sample ID (\"\") in `SUBJECT_SAMPLE_FACTOR` block.")
    assert type(validation_errors[0]) == type(test_error) and validation_errors[0].args == test_error.args
    test_error = ValueError("Sample with no Sample ID (\"\") in `MS_METABOLITE_DATA` block (usually caused by "
                            "extraneous TAB at the end of line).")
    assert type(validation_errors[1]) == type(test_error) and validation_errors[1].args == test_error.args
    test_error = ValueError("`MS_METABOLITE_DATA` block contains additional samples not found in "
                            "`SUBJECT_SAMPLE_FACTORS` block.\n\tAdditional samples: {'LabF_115873'}")
    assert type(validation_errors[2]) == type(test_error) and validation_errors[2].args == test_error.args


@pytest.mark.parametrize("file_source", [[
    "tests/example_data/validation_files/AN000041_error_1.txt",
    "tests/example_data/validation_files/AN000041_error_2.txt"
]])
def test_validate_nmr_samples(file_source):
    mwfile = next(mwtab.read_files(file_source[0]))
    validation_errors = mwtab.validate_file(mwfile, validate_factors=False, validate_features=False,
                                            validate_schema=False, validate_data=False, test=True)
    assert len(validation_errors) == 1
    test_error = KeyError("Missing key `Bin range(ppm)` in `NMR_BINNED_DATA` block.")
    assert type(validation_errors[0]) == type(test_error) and validation_errors[0].args == test_error.args

    mwfile = next(mwtab.read_files(file_source[1]))
    validation_errors = mwtab.validate_file(mwfile, validate_factors=False, validate_features=False,
                                            validate_schema=False, validate_data=False, test=True)
    assert len(validation_errors) == 3
    test_error = ValueError("Sample with no Sample ID (\"\") in `SUBJECT_SAMPLE_FACTOR` block.")
    assert type(validation_errors[0]) == type(test_error) and validation_errors[0].args == test_error.args
    test_error = ValueError("Sample with no sample ID (\"\") in `NMR_BINNED_DATA` block (usually caused by extraneous "
                            "TAB at the end of line).")
    assert type(validation_errors[1]) == type(test_error) and validation_errors[1].args == test_error.args
    test_error = ValueError("`NMR_BINNED_DATA` block contains additional samples not found in `SUBJECT_SAMPLE_FACTORS` "
                            "block.\n\tAdditional samples: {'C0559'}")
    assert type(validation_errors[2]) == type(test_error) and validation_errors[2].args == test_error.args


@pytest.mark.parametrize("file_source", [
    "tests/example_data/validation_files/AN000001_error_3.txt",

])
def test_validate_factors(file_source):
    mwfile = next(mwtab.read_files(file_source))
    validation_errors = mwtab.validate_file(mwfile, validate_samples=False, validate_features=False,
                                            validate_schema=False, validate_data=False, test=True)
    assert len(validation_errors) == 1
    test_error = KeyError("Missing key `Factors` in `MS_METABOLITE_DATA` block.")
    assert type(validation_errors[0]) == type(test_error) and validation_errors[0].args == test_error.args


@pytest.mark.parametrize("file_source", [[
    "tests/example_data/validation_files/AN000001_error_5.txt",
    "tests/example_data/validation_files/AN000001_error_6.txt"
]])
def test_validate_metabolites(file_source):
    mwfile = next(mwtab.read_files(file_source[0]))
    validation_errors = mwtab.validate_file(mwfile, validate_samples=False, validate_factors=False,
                                            validate_schema=False, validate_data=False, test=True)
    assert len(validation_errors) == 1
    test_error = KeyError("Missing key `metabolite_name` in `METABOLITES` block.")
    assert type(validation_errors[0]) == type(test_error) and validation_errors[0].args == test_error.args

    mwfile = next(mwtab.read_files(file_source[1]))
    validation_errors = mwtab.validate_file(mwfile, validate_samples=False, validate_factors=False,
                                            validate_schema=False, validate_data=False, test=True)
    assert len(validation_errors) == 4
    test_error = ValueError("Feature with no name (\"\") in `MS_METABOLITE_DATA` block.")
    assert type(validation_errors[0]) == type(test_error) and validation_errors[0].args == test_error.args
    test_error = ValueError("Feature with no name (\"\") in `METABOLITES` block.")
    assert type(validation_errors[1]) == type(test_error) and validation_errors[1].args == test_error.args
    test_error = ValueError("`MS_METABOLITE_DATA` block contains additional features not found in `METABOLITES` "
                            "block.\n\tAdditional features: ['1-monostearin', 'xylose']")
    assert type(validation_errors[2]) == type(test_error) and validation_errors[2].args == test_error.args
    test_error = ValueError("`METABOLITES` block contains additional features not found in `MS_METABOLITE_DATA` "
                            "block.\n\tAdditional features: ['1,2,4-benzenetriol', 'xylonic acid']")
    assert type(validation_errors[3]) == type(test_error) and validation_errors[3].args == test_error.args
