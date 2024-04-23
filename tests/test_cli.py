from __future__ import unicode_literals
import csv
import json
import mwtab
import os
import pytest
import shutil
import time
import pathlib
import subprocess


def teardown_module(module):
    path = pathlib.Path("tests/example_data/tmp/")
    if os.path.exists(path):
        shutil.rmtree(path)
        time_to_wait=10
        time_counter = 0
        while path.exists():
            time.sleep(1)
            time_counter += 1
            if time_counter > time_to_wait:
                raise FileExistsError(path + " was not deleted within " + str(time_to_wait) + " seconds, so it is assumed that it won't be and something went wrong.")


@pytest.mark.parametrize("files_source", [
    "204",
    "AN000204",
    "https://www.metabolomicsworkbench.org/rest/study/analysis_id/AN000204/mwtab/txt",
    "tests/example_data/mwtab_files/ST000122_AN000204.txt",
    "tests/example_data/mwtab_files/ST000122_AN000204.json",
    "tests/example_data/mwtab_files",
    "tests/example_data/mwtab_files.zip",
    "tests/example_data/mwtab_files.tar",
    "tests/example_data/mwtab_files.tar.gz",
    "tests/example_data/mwtab_files.tar.bz2"
])
def test_validate_command(files_source):
    command = "python -m mwtab validate {}".format(files_source)
    assert os.system(command) == 0


def test_validate_command_error_recovery():
    """Test that the validate command can have an error and move to the next one."""
    command = "python -m mwtab validate tests/example_data/files_to_test_error_recovery"
    command = command.split(" ")
    subp = subprocess.run(command, capture_output=True, encoding="UTF-8")
    assert subp.returncode == 0
    assert "1_no_error.json" in subp.stdout
    assert "2_error.json" in subp.stdout
    assert "3_no_error.txt" in subp.stdout
    assert "ValueError: Blank input string retrieved from source." in subp.stdout


@pytest.mark.parametrize("from_path, to_path, from_format, to_format", [
    # one-to-one file conversions
    ("tests/example_data/mwtab_files/ST000122_AN000204.txt", "tests/example_data/tmp/json/ST000122_AN000204.json", "mwtab", "json"),
    ("tests/example_data/mwtab_files/ST000122_AN000204.txt", "tests/example_data/tmp/json/ST000122_AN000204.json.gz", "mwtab", "json"),
    ("tests/example_data/mwtab_files/ST000122_AN000204.txt", "tests/example_data/tmp/json/ST000122_AN000204.json.bz2", "mwtab", "json"),
    ("tests/example_data/tmp/json/ST000122_AN000204.json", "tests/example_data/tmp/mwtab/ST000122_AN000204.txt", "json", "mwtab"),
    ("tests/example_data/tmp/json/ST000122_AN000204.json", "tests/example_data/tmp/mwtab/ST000122_AN000204.txt.gz", "json", "mwtab"),
    ("tests/example_data/tmp/json/ST000122_AN000204.json", "tests/example_data/tmp/mwtab/ST000122_AN000204.txt.bz2", "json", "mwtab"),
    ("tests/example_data/tmp/json/ST000122_AN000204.json.gz", "tests/example_data/tmp/mwtab/ST000122_AN000204.txt", "json", "mwtab"),
    ("tests/example_data/tmp/json/ST000122_AN000204.json.gz", "tests/example_data/tmp/mwtab/ST000122_AN000204.txt.gz", "json", "mwtab"),
    ("tests/example_data/tmp/json/ST000122_AN000204.json.gz", "tests/example_data/tmp/mwtab/ST000122_AN000204.txt.bz2", "json", "mwtab"),
    ("tests/example_data/tmp/json/ST000122_AN000204.json.bz2", "tests/example_data/tmp/mwtab/ST000122_AN000204.txt", "json", "mwtab"),
    ("tests/example_data/tmp/json/ST000122_AN000204.json.bz2", "tests/example_data/tmp/mwtab/ST000122_AN000204.txt.gz", "json", "mwtab"),
    ("tests/example_data/tmp/json/ST000122_AN000204.json.bz2", "tests/example_data/tmp/mwtab/ST000122_AN000204.txt.bz2", "json", "mwtab"),
    ("tests/example_data/mwtab_files/ST000122_AN000204.json", "tests/example_data/tmp/mwtab/ST000122_AN000204.txt", "json", "mwtab"),
    ("tests/example_data/mwtab_files/ST000122_AN000204.json", "tests/example_data/tmp/mwtab/ST000122_AN000204.txt.gz", "json", "mwtab"),
    ("tests/example_data/mwtab_files/ST000122_AN000204.json", "tests/example_data/tmp/mwtab/ST000122_AN000204.txt.bz2", "json", "mwtab"),
    ("tests/example_data/tmp/mwtab/ST000122_AN000204.txt", "tests/example_data/tmp/json/ST000122_AN000204.json", "mwtab", "json"),
    ("tests/example_data/tmp/mwtab/ST000122_AN000204.txt", "tests/example_data/tmp/json/ST000122_AN000204.json.gz", "mwtab", "json"),
    ("tests/example_data/tmp/mwtab/ST000122_AN000204.txt", "tests/example_data/tmp/json/ST000122_AN000204.json.bz2", "mwtab", "json"),
    ("tests/example_data/tmp/mwtab/ST000122_AN000204.txt.gz", "tests/example_data/tmp/json/ST000122_AN000204.json", "mwtab", "json"),
    ("tests/example_data/tmp/mwtab/ST000122_AN000204.txt.gz", "tests/example_data/tmp/json/ST000122_AN000204.json.gz", "mwtab", "json"),
    ("tests/example_data/tmp/mwtab/ST000122_AN000204.txt.gz", "tests/example_data/tmp/json/ST000122_AN000204.json.bz2", "mwtab", "json"),
    ("tests/example_data/tmp/mwtab/ST000122_AN000204.txt.bz2", "tests/example_data/tmp/json/ST000122_AN000204.json", "mwtab", "json"),
    ("tests/example_data/tmp/mwtab/ST000122_AN000204.txt.bz2", "tests/example_data/tmp/json/ST000122_AN000204.json.gz", "mwtab", "json"),
    ("tests/example_data/tmp/mwtab/ST000122_AN000204.txt.bz2", "tests/example_data/tmp/json/ST000122_AN000204.json.bz2", "mwtab", "json"),
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
    mwtabfiles_list = [file for file in mwtabfile_generator]
    mwtabfiles_study_ids_set = set(mwf.study_id for mwf in mwtabfiles_list)
    mwtabfiles_analysis_ids_set = set(mwf.analysis_id for mwf in mwtabfiles_list)
    assert mwtabfiles_study_ids_set.issubset({"ST000122"})
    assert mwtabfiles_analysis_ids_set.issubset({"AN000204"})


def test_convert_command_mwrest_option():
    """Test that the --mw-rest option works for the convert command. Test this by setting a bad url and making sure there is an error."""
    command = "python -m mwtab convert 204 tests/example_data/tmp/tmp.json --mw-rest=https://bad_url --from-format=mwtab --to-format=json"
    command = command.split(" ")
    subp = subprocess.run(command, capture_output=True, encoding="UTF-8")
    assert "urlopen error" in subp.stdout
    assert subp.returncode == 0


def test_convert_command_error_recovery():
    """Test that the convert command can recover from a bad file."""
    command = "python -m mwtab convert tests/example_data/files_to_test_error_recovery tests/example_data/tmp/ --from-format=mwtab --to-format=json"
    command = command.split(" ")
    subp = subprocess.run(command, capture_output=True, encoding="UTF-8")
    assert subp.returncode == 0
    assert "1_no_error.json" in subp.stdout
    assert "2_error.json" in subp.stdout
    assert "3_no_error.txt" in subp.stdout
    assert "ValueError: Blank input string retrieved from source." in subp.stdout


@pytest.mark.parametrize("command", [
    # download by url
    "python -m mwtab download url https://www.metabolomicsworkbench.org/rest/study/study_id/ST000001/summary --to-path=tests/example_data/tmp/tmp.txt",
    # download by study methods
    "python -m mwtab download study 2 --to-path=tests/example_data/tmp/tmp.txt --output-format=txt",
    "python -m mwtab download study ST000002 --to-path=tests/example_data/tmp/tmp.txt --output-format=txt",
    "python -m mwtab download study study_id ST000002 summary --to-path=tests/example_data/tmp/tmp.txt",
    "python -m mwtab download study study_id ST analysis --to-path=tests/example_data/tmp/tmp.txt",
    # download compound | refmet | gene | protein
    "python -m mwtab download compound regno 11 name --to-path=tests/example_data/tmp/tmp.txt",
    "python -m mwtab download refmet name Cholesterol all --to-path=tests/example_data/tmp/tmp.txt",
    "python -m mwtab download gene gene_symbol acaca all --to-path=tests/example_data/tmp/tmp.txt",
    "python -m mwtab download protein uniprot_id Q13085 all --to-path=tests/example_data/tmp/tmp.txt",
    # download moverz
    "python -m mwtab download moverz MB 635.52 M+H 0.5 --to-path=tests/example_data/tmp/tmp.txt",
    "python -m mwtab download moverz LIPIDS 513.45 M-2H 0.2 --to-path=tests/example_data/tmp/tmp.txt",
    "python -m mwtab download moverz REFMET 255.2 M+H 0.2 --to-path=tests/example_data/tmp/tmp.txt",
    # download exactmass
    "python -m mwtab download exactmass \"PC(34:1)\" M+H --to-path=tests/example_data/tmp/tmp.txt",
    "python -m mwtab download exactmass  \"GlcCer(d42:2)\" M-H --to-path=tests/example_data/tmp/tmp.txt",

])
def test_download_command(command):
    assert os.system(command) == 0

    file_str = ""
    with open("tests/example_data/tmp/tmp.txt", "r") as fh:
        file_str = fh.read()
        fh.close()
    with open("tests/example_data/tmp/tmp.txt", "w") as fh:
        fh.close()
    assert file_str


def test_download_url_command_error():
    """Test that the download url command can have an error and print appropriately."""
    command = "python -m mwtab download url https://bad_url --to-path=tests/example_data/tmp/tmp.txt"
    command = command.split(" ")
    subp = subprocess.run(command, capture_output=True, encoding="UTF-8")
    assert subp.returncode == 1


def test_download_study_error_recovery():
    """Test that the download study command can recover from an error."""
    command = "python -m mwtab download study tests/example_data/other_test_data/download_studies_with_error.json --to-path=tests/example_data/tmp/"
    command = command.split(" ")
    subp = subprocess.run(command, capture_output=True, encoding="UTF-8")
    assert "AN9999999 could not be downloaded" in subp.stdout
    assert ('When trying to download a file for the value, "AN000206", '
            'a blank file or an error was returned, so no file was created for it.') in subp.stdout
    assert "ValueError: Invalid Metabolomics Workbench analysis ID for a study (AN<6-digit integer>)" in subp.stdout
    assert subp.returncode == 0


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
                assert set(data[0]) == {"SUBJECT_TYPE", "Human"}
            else:
                assert set(data[0]) == {"metadata", "value0"}
                assert set(data[1]) == {"SUBJECT_TYPE", "Human"}
        elif to_format == "json":
            data = json.load(f)
            data["SUBJECT_TYPE"] = set(data["SUBJECT_TYPE"])
            assert data == {"SUBJECT_TYPE": {"Human"}}
        else:
            assert False


def test_extract_metadata_error_recovery():
    """Test that the extract metadata command can recover from an error."""
    command = "python -m mwtab extract metadata tests/example_data/files_to_test_error_recovery tests/example_data/tmp/ SUBJECT_TYPE"
    command = command.split(" ")
    subp = subprocess.run(command, capture_output=True, encoding="UTF-8")
    assert subp.returncode == 0
    assert "1_no_error.json" in subp.stdout
    assert "2_error.json" in subp.stdout
    assert "3_no_error.txt" in subp.stdout
    assert "ValueError: Blank input string retrieved from source." in subp.stdout


def test_extract_metadata_no_data_extracted():
    """Test that the extract metadata command prints a message when no data is extracted and no file is made."""
    command = "python -m mwtab extract metadata tests/example_data/mwtab_files/ tests/example_data/tmp/ asdf"
    command = command.split(" ")
    subp = subprocess.run(command, capture_output=True, encoding="UTF-8")
    assert not os.path.exists("tests/example_data/tmp/")
    assert subp.returncode == 0
    assert "No metadata extracted. No file was saved." in subp.stdout


@pytest.mark.parametrize("from_path, to_path, key, value, to_format, no_header", [
    ("tests/example_data/mwtab_files/", "tests/example_data/tmp/test_extract_metabolites", "SU:SUBJECT_TYPE", "Human", "csv", " --no-header"),
    ("tests/example_data/mwtab_files/", "tests/example_data/tmp/test_extract_metabolites", "SU:SUBJECT_TYPE", "Human", "csv", ""),
    ("tests/example_data/mwtab_files/", "tests/example_data/tmp/test_extract_metabolites", "SU:SUBJECT_TYPE", "Human", "json", ""),
    ("tests/example_data/mwtab_files/", "tests/example_data/tmp/test_extract_metabolites.csv", "SU:SUBJECT_TYPE", "Human", "csv", ""),
    ("tests/example_data/mwtab_files/", "tests/example_data/tmp/test_extract_metabolites.json", "SU:SUBJECT_TYPE", "Human", "json", ""),
    ("tests/example_data/mwtab_files/", "tests/example_data/tmp/test_extract_metabolites", "SU:SUBJECT_TYPE", "\"r'(Human)'\"", "csv", " --no-header"),
    ("tests/example_data/mwtab_files/", "tests/example_data/tmp/test_extract_metabolites", "SU:SUBJECT_TYPE", "\"r'(Human)'\"", "csv", ""),
    ("tests/example_data/mwtab_files/", "tests/example_data/tmp/test_extract_metabolites", "SU:SUBJECT_TYPE", "\"r'(Human)'\"", "json", "")
])
def test_extract_metabolites_command(from_path, to_path, key, value, to_format, no_header):
    command = "python -m mwtab extract metabolites {} {} {} {} --to-format={}{}".format(
        from_path, to_path, key, value, to_format, no_header
    )
    assert os.system(command) == 0
    
    expected_data = [['17-hydroxypregnenolone', '1', '1', '13'],
                     ['17-hydroxyprogesterone', '1', '1', '29'],
                     ['Allodihydrotestosterone', '1', '1', '18'],
                     ['Androstenedione', '1', '1', '42'],
                     ['Androstenolone (DHEA)', '1', '1', '42'],
                     ['Cortexolone', '1', '1', '21'],
                     ['Cortexone', '1', '1', '41'],
                     ['Corticosterone_ DOC', '1', '1', '29'],
                     ['Cortisol', '1', '1', '42'],
                     ['Estradiol', '1', '1', '42'],
                     ['Estrone', '1', '1', '42'],
                     ['Pregnenolone', '1', '1', '12'],
                     ['Progesterone', '1', '1', '39'],
                     ['Testosterone', '1', '1', '42']]

    if to_format == "csv":
        filepath = to_path
        if not os.path.splitext(filepath)[1]:
            filepath += ".csv"
        with open(filepath, "r") as fh:
            data = list(csv.reader(fh))
            if bool(no_header):
                assert data == expected_data
            else:
                assert data[0] == ["metabolite_name", "num-studies", "num_analyses", "num_samples"]
                assert data[1:] == expected_data
            fh.close()
    elif to_format == 'json':
        filepath = to_path
        if not os.path.splitext(filepath)[1]:
            filepath += ".json"
        with open(filepath, "r") as fh:
            text = fh.read()
            fh.close()
        assert text
    else:
        assert False


def test_extract_metabolites_error_recovery():
    """Test that the extract metabolites command can recover from an error."""
    command = "python -m mwtab extract metabolites tests/example_data/files_to_test_error_recovery tests/example_data/tmp/ SU:SUBJECT_TYPE Human"
    command = command.split(" ")
    subp = subprocess.run(command, capture_output=True, encoding="UTF-8")
    # print()
    # print("Output:")
    # print(subp.stdout)
    # print()
    # print("Error:")
    # print(subp.stderr)
    assert subp.returncode == 0
    assert "1_no_error.json" in subp.stdout
    assert "2_error.json" in subp.stdout
    assert "3_no_error.txt" in subp.stdout
    assert "ValueError: Blank input string retrieved from source." in subp.stdout


def test_extract_metabolites_no_data_extracted():
    """Test that the extract metabolites command prints a message when no data is extracted and no file is made."""
    command = "python -m mwtab extract metabolites tests/example_data/mwtab_files/ tests/example_data/tmp/ SU:SUBJECT_TYPE Plant"
    command = command.split(" ")
    subp = subprocess.run(command, capture_output=True, encoding="UTF-8")
    assert not os.path.exists("tests/example_data/tmp/")
    assert subp.returncode == 0
    assert ("No metabolites extracted. No file was saved. " 
            "This is likely due to key value pairs filtering all of the studies out. "
            "Check your key value pairs and data.") in subp.stdout


