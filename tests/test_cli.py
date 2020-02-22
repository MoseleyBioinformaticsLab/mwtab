from __future__ import unicode_literals
import csv
import json
import mwtab
import os
import pytest
import shutil


def teardown_module(module):
    if os.path.exists("tests/example_data/tmp/"):
        shutil.rmtree("tests/example_data/tmp")


@pytest.mark.parametrize("files_source", [
    "2",  # changed to 2 from 1. As of 02/12/2020 AN000001 contains an error preventing it from being parsed.
    "tests/example_data/mwtab_files",
    "tests/example_data/mwtab_files.zip",
    "tests/example_data/mwtab_files.tar.gz",
    "tests/example_data/mwtab_files.tar.bz2"
])
def test_validate_command(files_source):
    command = "python -m mwtab validate {}".format(files_source)
    assert os.system(command) == 0


@pytest.mark.parametrize("from_path, to_path, from_format, to_format", [
    # one-to-one file conversions
    ("tests/example_data/mwtab_files/ST000001_AN000001.txt", "tests/example_data/tmp/json/ST000001_AN000001.json", "mwtab", "json"),
    ("tests/example_data/mwtab_files/ST000001_AN000001.txt", "tests/example_data/tmp/json/ST000001_AN000001.json.gz", "mwtab", "json"),
    ("tests/example_data/mwtab_files/ST000001_AN000001.txt", "tests/example_data/tmp/json/ST000001_AN000001.json.bz2", "mwtab", "json"),
    ("tests/example_data/tmp/json/ST000001_AN000001.json.gz", "tests/example_data/tmp/mwtab/ST000001_AN000001.txt", "json", "mwtab"),
    ("tests/example_data/tmp/json/ST000001_AN000001.json.gz", "tests/example_data/tmp/mwtab/ST000001_AN000001.txt.gz", "json", "mwtab"),
    ("tests/example_data/tmp/json/ST000001_AN000001.json.gz", "tests/example_data/tmp/mwtab/ST000001_AN000001.txt.bz2", "json", "mwtab"),
    ("tests/example_data/tmp/json/ST000001_AN000001.json.bz2", "tests/example_data/tmp/mwtab/ST000001_AN000001.txt", "json", "mwtab"),
    ("tests/example_data/tmp/json/ST000001_AN000001.json.bz2", "tests/example_data/tmp/mwtab/ST000001_AN000001.txt.gz", "json", "mwtab"),
    ("tests/example_data/tmp/json/ST000001_AN000001.json.bz2", "tests/example_data/tmp/mwtab/ST000001_AN000001.txt.bz2", "json", "mwtab"),
    # many-to-many file conversions
    ("tests/example_data/mwtab_files", "tests/example_data/tmp/json/dir/mwtab_files_json", "mwtab", "json"),
    ("tests/example_data/mwtab_files", "tests/example_data/tmp/json/dir/mwtab_files_json.zip", "mwtab", "json"),
    ("tests/example_data/mwtab_files", "tests/example_data/tmp/json/dir/mwtab_files_json.tar", "mwtab", "json"),
    ("tests/example_data/mwtab_files", "tests/example_data/tmp/json/dir/mwtab_files_json.tar.gz", "mwtab", "json"),
    ("tests/example_data/mwtab_files", "tests/example_data/tmp/json/dir/mwtab_files_json.tar.bz2", "mwtab", "json"),
    ("tests/example_data/tmp/json/dir/mwtab_files_json.zip", "tests/example_data/tmp/mwtab/zip/mwtab_files_mwtab", "json", "mwtab"),
    ("tests/example_data/tmp/json/dir/mwtab_files_json.zip", "tests/example_data/tmp/mwtab/zip/mwtab_files_mwtab.zip", "json", "mwtab"),
    ("tests/example_data/tmp/json/dir/mwtab_files_json.zip", "tests/example_data/tmp/mwtab/zip/mwtab_files_mwtab.tar", "json", "mwtab"),
    ("tests/example_data/tmp/json/dir/mwtab_files_json.zip", "tests/example_data/tmp/mwtab/zip/mwtab_files_mwtab.tar.gz", "json", "mwtab"),
    ("tests/example_data/tmp/json/dir/mwtab_files_json.zip", "tests/example_data/tmp/mwtab/zip/mwtab_files_mwtab.tar.bz2", "json", "mwtab"),
    ("tests/example_data/tmp/json/dir/mwtab_files_json.tar", "tests/example_data/tmp/mwtab/tar/mwtab_files_mwtab", "json", "mwtab"),
    ("tests/example_data/tmp/json/dir/mwtab_files_json.tar", "tests/example_data/tmp/mwtab/tar/mwtab_files_mwtab.zip", "json", "mwtab"),
    ("tests/example_data/tmp/json/dir/mwtab_files_json.tar", "tests/example_data/tmp/mwtab/tar/mwtab_files_mwtab.tar", "json", "mwtab"),
    ("tests/example_data/tmp/json/dir/mwtab_files_json.tar", "tests/example_data/tmp/mwtab/tar/mwtab_files_mwtab.tar.gz", "json", "mwtab"),
    ("tests/example_data/tmp/json/dir/mwtab_files_json.tar", "tests/example_data/tmp/mwtab/tar/mwtab_files_mwtab.tar.bz2", "json", "mwtab"),
    ("tests/example_data/tmp/json/dir/mwtab_files_json.tar.gz", "tests/example_data/tmp/mwtab/targz/mwtab_files_mwtab", "json", "mwtab"),
    ("tests/example_data/tmp/json/dir/mwtab_files_json.tar.gz", "tests/example_data/tmp/mwtab/targz/mwtab_files_mwtab.zip", "json", "mwtab"),
    ("tests/example_data/tmp/json/dir/mwtab_files_json.tar.gz", "tests/example_data/tmp/mwtab/targz/mwtab_files_mwtab.tar", "json", "mwtab"),
    ("tests/example_data/tmp/json/dir/mwtab_files_json.tar.gz", "tests/example_data/tmp/mwtab/targz/mwtab_files_mwtab.tar.gz", "json", "mwtab"),
    ("tests/example_data/tmp/json/dir/mwtab_files_json.tar.gz", "tests/example_data/tmp/mwtab/targz/mwtab_files_mwtab.tar.bz2", "json", "mwtab"),
    ("tests/example_data/tmp/json/dir/mwtab_files_json.tar.bz2", "tests/example_data/tmp/mwtab/tarbz2/mwtab_files_mwtab", "json", "mwtab"),
    ("tests/example_data/tmp/json/dir/mwtab_files_json.tar.bz2", "tests/example_data/tmp/mwtab/tarbz2/mwtab_files_mwtab.zip", "json", "mwtab"),
    ("tests/example_data/tmp/json/dir/mwtab_files_json.tar.bz2", "tests/example_data/tmp/mwtab/tarbz2/mwtab_files_mwtab.tar", "json", "mwtab"),
    ("tests/example_data/tmp/json/dir/mwtab_files_json.tar.bz2", "tests/example_data/tmp/mwtab/tarbz2/mwtab_files_mwtab.tar.gz", "json", "mwtab"),
    ("tests/example_data/tmp/json/dir/mwtab_files_json.tar.bz2", "tests/example_data/tmp/mwtab/tarbz2/mwtab_files_mwtab.tar.bz2", "json", "mwtab")
])
def test_convert_command(from_path, to_path, from_format, to_format):
    command = "python -m mwtab convert {} {} --from-format={} --to-format={}".format(
        from_path, to_path, from_format, to_format
    )
    assert os.system(command) == 0

    mwtabfile_generator = mwtab.read_files(to_path)
    mwtabfiles_list = list(mwtabfile_generator)
    mwtabfiles_study_ids_set = set(mwf.study_id for mwf in mwtabfiles_list)
    mwtabfiles_analysis_ids_set = set(mwf.analysis_id for mwf in mwtabfiles_list)
    assert mwtabfiles_study_ids_set.issubset({"ST000001", "ST000002"})
    assert mwtabfiles_analysis_ids_set.issubset({"AN000001", "AN000002"})


@pytest.mark.parametrize("input_value, to_path", [
    ("AN000002", "tests/example_data/tmp")
])
def test_download_command(input_value, to_path):
    command = "python -m mwtab download {} --to-path={}".format(input_value, to_path)
    assert os.system(command) == 0

    mwtabfile = next(mwtab.read_files(to_path))
    assert mwtabfile.study_id == "ST000002"
    assert mwtabfile.analysis_id == "AN000002"


@pytest.mark.parametrize("from_path, to_path, key, to_format, no_header", [
    ("tests/example_data/mwtab_files/", "tests/example_data/tmp/test_extract_metadata", "SUBJECT_TYPE", "csv", " --no-header"),
    ("tests/example_data/mwtab_files/", "tests/example_data/tmp/test_extract_metadata", "SUBJECT_TYPE", "csv", ""),
    ("tests/example_data/mwtab_files/", "tests/example_data/tmp/test_extract_metadata", "SUBJECT_TYPE", "json", "")
])
def test_extract_metadata_command(from_path, to_path, key, to_format, no_header):
    command = "python -m mwtab extract metadata {} {} {} --to-format={}{}".format(
        from_path, to_path, key, to_format, no_header
    )
    assert os.system(command) == 0

    with open(".".join([to_path, to_format]), "r") as f:
        if to_format == "csv":
            data = list(csv.reader(f))
            if bool(no_header):
                assert set(data[0]) == {"SUBJECT_TYPE", "Human", "Plant"}
            else:
                assert set(data[0]) == {"metadata", "value0", "value1"}
                assert set(data[1]) == {"SUBJECT_TYPE", "Human", "Plant"}
        elif to_format == "json":
            data = json.load(f)
            data["SUBJECT_TYPE"] = set(data["SUBJECT_TYPE"])
            assert data == {"SUBJECT_TYPE": {"Human", "Plant"}}
        else:
            assert False


@pytest.mark.parametrize("from_path, to_path, key, value, to_format, no_header", [
    ("tests/example_data/mwtab_files/", "tests/example_data/tmp/test_extract_metabolites", "SU:SUBJECT_TYPE", "Plant", "csv", " --no-header"),
    ("tests/example_data/mwtab_files/", "tests/example_data/tmp/test_extract_metabolites", "SU:SUBJECT_TYPE", "Plant", "csv", ""),
    ("tests/example_data/mwtab_files/", "tests/example_data/tmp/test_extract_metabolites", "SU:SUBJECT_TYPE", "Plant", "json", ""),
    ("tests/example_data/mwtab_files/", "tests/example_data/tmp/test_extract_metabolites", "SU:SUBJECT_TYPE", "\"r'(Plant)'\"", "csv", " --no-header"),
    ("tests/example_data/mwtab_files/", "tests/example_data/tmp/test_extract_metabolites", "SU:SUBJECT_TYPE", "\"r'(Plant)'\"", "csv", ""),
    ("tests/example_data/mwtab_files/", "tests/example_data/tmp/test_extract_metabolites", "SU:SUBJECT_TYPE", "\"r'(Plant)'\"", "json", "")
])
def test_extract_metabolites_command(from_path, to_path, key, value, to_format, no_header):
    command = "python -m mwtab extract metabolites {} {} {} {} --to-format={}{}".format(
        from_path, to_path, key, value, to_format, no_header
    )
    print(command)
    assert os.system(command) == 0

    if to_format == "csv":
        with open(".".join([to_path, to_format]), "r") as f:
            data = list(csv.reader(f))
            if bool(no_header):
                assert set(data[0]) == {"1,2,4-benzenetriol", "1", "1", "24"}
                assert len(data) == 191
            else:
                assert set(data[0]) == {"metabolite_name", "num-studies", "num_analyses", "num_samples"}
                # assert set(data[1]) == {"1,2,4-benzenetriol", "1", "1", "24"}
                assert len(data) == 192
