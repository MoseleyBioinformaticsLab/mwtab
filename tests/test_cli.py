from __future__ import unicode_literals
import csv
import json
import mwtab
import os
import pytest
import time
import pathlib
import subprocess

from fixtures import teardown_module

from mwtab import cli

teardown_module = teardown_module

def delete_ST000001_cwd():
    path = pathlib.Path(os.path.join(os.getcwd(), 'https%3A%2F%2Fwww_metabolomicsworkbench_org%2Frest%2Fstudy%2Fstudy_id%2FST000001%2Fsummary.txt'))
    if path.exists():
        os.remove(path)
        time_to_wait=10
        time_counter = 0
        while path.exists():
            time.sleep(1)
            time_counter += 1
            if time_counter > time_to_wait:
                raise FileExistsError(path + " was not deleted within " + str(time_to_wait) + " seconds, so it is assumed that it won't be and something went wrong.")

@pytest.fixture()
def teardown_module_cwd():
    delete_ST000001_cwd()
    yield
    delete_ST000001_cwd()
    

@pytest.fixture()
def disable_network_calls(monkeypatch):
    def stunted_get():
        raise RuntimeError("Network access not allowed during testing!")
    monkeypatch.setattr("urllib.request.urlopen", lambda *args, **kwargs: stunted_get())
    monkeypatch.setattr("urllib.parse.urlparse", lambda *args, **kwargs: stunted_get())
    monkeypatch.setattr("mwtab.fileio.urlopen", lambda *args, **kwargs: stunted_get())
    monkeypatch.setattr("mwtab.fileio.urlparse", lambda *args, **kwargs: stunted_get())

@pytest.fixture()
def disable_sleep(monkeypatch):
    def no_sleep(arg):
        pass
    monkeypatch.setattr('mwtab.cli.time.sleep', no_sleep)


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

@pytest.mark.parametrize("file_sources", [
    ("tests/example_data/mwtab_files/ST000122_AN000204.txt", 'tests/example_data/tmp/validation_errors.json'),
    ("tests/example_data/mwtab_files/ST000122_AN000204.txt", 'tests/example_data/tmp')
])
def test_validate_command_saving(file_sources, teardown_module):
    command = "python -m mwtab validate {} --to-path={}".format(*file_sources)
    assert os.system(command) == 0
    assert os.path.exists(file_sources[1])


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


def test_convert_command_mwrest_option(teardown_module):
    """Test that the --mw-rest option works for the convert command. Test this by setting a bad url and making sure there is an error."""
    command = "python -m mwtab convert 204 tests/example_data/tmp/tmp.json --mw-rest=https://bad_url --from-format=mwtab --to-format=json"
    command = command.split(" ")
    subp = subprocess.run(command, capture_output=True, encoding="UTF-8")
    assert "urlopen error" in subp.stdout
    assert subp.returncode == 0


def test_convert_command_error_recovery(teardown_module):
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
def test_download_command(command, teardown_module):
    assert os.system(command) == 0

    file_str = ""
    with open("tests/example_data/tmp/tmp.txt", "r") as fh:
        file_str = fh.read()
        fh.close()
    with open("tests/example_data/tmp/tmp.txt", "w") as fh:
        fh.close()
    assert file_str


def test_download_command_cwd(teardown_module_cwd):
    assert os.system('python -m mwtab download url https://www.metabolomicsworkbench.org/rest/study/study_id/ST000001/summary') == 0
    file_path = os.path.join(os.getcwd(), 'https%3A%2F%2Fwww_metabolomicsworkbench_org%2Frest%2Fstudy%2Fstudy_id%2FST000001%2Fsummary.txt')

    file_str = ""
    with open(file_path, "r") as fh:
        file_str = fh.read()
        fh.close()
    with open(file_path, "w") as fh:
        fh.close()
    assert file_str


def test_download_all_analysis_ids(teardown_module, disable_network_calls, disable_sleep, capsys, mocker):
    cmdargs = {
        '--verbose': True,
        '--force': False,
        '--silent': False,
        '--mw-rest': 'https://www.metabolomicsworkbench.org/rest/',
        'convert': False,
        'validate': False,
        'download': True,
        '<url>': None,
        'study': True,
        'all': True,
        '--input-item': 'analysis_id'
        }
    mocker.patch('mwtab.cli.mwrest.analysis_ids', side_effect = [['AN000001', 'AN000002']])
    mocker.patch('mwtab.cli.download_and_save_mwrest_file', side_effect = [None, None])
    
    cli.cli(cmdargs)
    captured = capsys.readouterr()
    assert 'AN000001' in captured.out
    assert 'AN000002' in captured.out

def test_download_all_analysis_ids2(teardown_module, disable_network_calls, disable_sleep, capsys, mocker):
    """The same as the previous test, but with --input-item being null to test the default."""
    cmdargs = {
        '--verbose': True,
        '--force': False,
        '--silent': False,
        '--mw-rest': 'https://www.metabolomicsworkbench.org/rest/',
        'convert': False,
        'validate': False,
        'download': True,
        '<url>': None,
        'study': True,
        'all': True,
        '--input-item': 'analysis_id'
        }
    mocker.patch('mwtab.cli.mwrest.analysis_ids', side_effect = [['AN000001', 'AN000002']])
    mocker.patch('mwtab.cli.download_and_save_mwrest_file', side_effect = [None, None])
    
    cli.cli(cmdargs)
    captured = capsys.readouterr()
    assert 'AN000001' in captured.out
    assert 'AN000002' in captured.out

def test_download_all_study_ids(teardown_module, disable_network_calls, disable_sleep, capsys, mocker):
    cmdargs = {
        '--verbose': True,
        '--force': False,
        '--silent': False,
        '--mw-rest': 'https://www.metabolomicsworkbench.org/rest/',
        'convert': False,
        'validate': False,
        'download': True,
        '<url>': None,
        'study': True,
        'all': True,
        '--input-item': 'study_id'
        }
    mocker.patch('mwtab.cli.mwrest.study_ids', side_effect = [['ST000001', 'ST000002']])
    mocker.patch('mwtab.cli.download_and_save_mwrest_file', side_effect = [None, None])
    
    cli.cli(cmdargs)
    captured = capsys.readouterr()
    assert 'ST000001' in captured.out
    assert 'ST000002' in captured.out

def test_download_all_download_error(teardown_module, disable_network_calls, disable_sleep, capsys, mocker):
    cmdargs = {
        '--verbose': True,
        '--force': False,
        '--silent': False,
        '--mw-rest': 'https://www.metabolomicsworkbench.org/rest/',
        'convert': False,
        'validate': False,
        'download': True,
        '<url>': None,
        'study': True,
        'all': True,
        '--input-item': 'analysis_id'
        }
    mocker.patch('mwtab.cli.mwrest.analysis_ids', side_effect = [['AN000001', 'AN000002']])
    mocker.patch('mwtab.cli.download_and_save_mwrest_file', side_effect = [Exception, None])
    
    cli.cli(cmdargs)
    captured = capsys.readouterr()
    assert 'AN000001' in captured.out
    assert 'AN000002' in captured.out
    assert 'Something went wrong and AN000001 could not be downloaded.' in captured.out

def test_download_all_bad_input_item(teardown_module, disable_network_calls, disable_sleep):
    cmdargs = {
        '--verbose': True,
        '--force': False,
        '--silent': False,
        '--mw-rest': 'https://www.metabolomicsworkbench.org/rest/',
        'convert': False,
        'validate': False,
        'download': True,
        '<url>': None,
        'study': True,
        'all': True,
        '--input-item': 'asdf'
        }
    
    with pytest.raises(ValueError, match = r'Unknown "--input-item" asdf'):
        cli.cli(cmdargs)

def test_download_study_file_list_analysis(teardown_module, disable_network_calls, disable_sleep, capsys, mocker):
    cmdargs = {
        '--verbose': True,
        '--force': False,
        '--silent': False,
        '--mw-rest': 'https://www.metabolomicsworkbench.org/rest/',
        'convert': False,
        'validate': False,
        'download': True,
        '<url>': None,
        'study': True,
        'all': None,
        '--input-item': 'analysis_id',
        '<input-item>': None,
        '<input-value>': 'tests/example_data/other_test_data/download_studies_an_ids.json'
        }
    mocker.patch('mwtab.cli.download_and_save_mwrest_file', side_effect = [None, None])
    
    cli.cli(cmdargs)
    captured = capsys.readouterr()
    assert 'Found 2 Files to be Downloaded' in captured.out
    assert 'AN000001' in captured.out
    assert 'AN000002' in captured.out

def test_download_study_file_list_study(teardown_module, disable_network_calls, disable_sleep, capsys, mocker):
    """Same as previous test, but the list is study IDs instead of AN IDs."""
    cmdargs = {
        '--verbose': True,
        '--force': False,
        '--silent': False,
        '--mw-rest': 'https://www.metabolomicsworkbench.org/rest/',
        'convert': False,
        'validate': False,
        'download': True,
        '<url>': None,
        'study': True,
        'all': None,
        '--input-item': 'study_id',
        '<input-item>': None,
        '<input-value>': 'tests/example_data/other_test_data/download_studies_study_ids.json'
        }
    mocker.patch('mwtab.cli.download_and_save_mwrest_file', side_effect = [None, None])
    
    cli.cli(cmdargs)
    captured = capsys.readouterr()
    assert 'Found 2 Files to be Downloaded' in captured.out
    assert 'ST000001' in captured.out
    assert 'ST000002' in captured.out


def test_download_study_file_list_mixed(teardown_module, disable_network_calls, disable_sleep, capsys, mocker):
    """Same as previous test, but the list values are mixed."""
    cmdargs = {
        '--verbose': True,
        '--force': False,
        '--silent': False,
        '--mw-rest': 'https://www.metabolomicsworkbench.org/rest/',
        'convert': False,
        'validate': False,
        'download': True,
        '<url>': None,
        'study': True,
        'all': None,
        '--input-item': None,
        '<input-item>': None,
        '<input-value>': 'tests/example_data/other_test_data/download_studies_mixed_ids.json'
        }
    mocker.patch('mwtab.cli.download_and_save_mwrest_file', side_effect = [None, None, None, Exception])
    
    cli.cli(cmdargs)
    captured = capsys.readouterr()
    assert 'Found 4 Files to be Downloaded' in captured.out
    assert '000001' in captured.out
    assert 'AN000002' in captured.out
    assert 'ST000003' in captured.out
    assert 'Something went wrong and asdf could not be downloaded.' in captured.out


def test_download_study_file_list_bad_input_item(teardown_module, disable_network_calls, disable_sleep):
    cmdargs = {
        '--verbose': True,
        '--force': False,
        '--silent': False,
        '--mw-rest': 'https://www.metabolomicsworkbench.org/rest/',
        'convert': False,
        'validate': False,
        'download': True,
        '<url>': None,
        'study': True,
        'all': None,
        '--input-item': 'asdf',
        '<input-item>': None,
        '<input-value>': 'tests/example_data/other_test_data/download_studies_mixed_ids.json'
        }
    
    with pytest.raises(ValueError, match = r'Unknown "--input-item" asdf'):
        cli.cli(cmdargs)


def test_download_study_file_list_bad_to_path(teardown_module, disable_network_calls, disable_sleep, capsys):
    cmdargs = {
        '--verbose': True,
        '--force': False,
        '--silent': False,
        '--mw-rest': 'https://www.metabolomicsworkbench.org/rest/',
        'convert': False,
        'validate': False,
        'download': True,
        '<url>': None,
        'study': True,
        'all': None,
        '--to-path': 'filename.txt',
        '<input-item>': None,
        '<input-value>': 'tests/example_data/other_test_data/download_studies_mixed_ids.json'
        }
    
    with pytest.raises(SystemExit):
        cli.cli(cmdargs)
    captured = capsys.readouterr()
    assert ('Error: The given "--to-path" option is a path to a file, '
            'but the given list of IDs to download is greater '
            'than 1. Please specify a directory to save the files to.') in captured.out


def test_download_url_command_error(teardown_module):
    """Test that the download url command can have an error and print appropriately."""
    command = "python -m mwtab download url https://bad_url --to-path=tests/example_data/tmp/tmp.txt"
    command = command.split(" ")
    subp = subprocess.run(command, capture_output=True, encoding="UTF-8")
    assert subp.returncode == 1


def test_download_study_error_recovery(teardown_module):
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
    ("tests/example_data/mwtab_files/", "tests/example_data/tmp/test_extract_metadata", "SUBJECT_TYPE", "json", ""),
])
def test_extract_metadata_command(from_path, to_path, key, to_format, no_header, teardown_module):
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


def test_extract_metadata_stdout():
    command = "python -m mwtab extract metadata tests/example_data/mwtab_files/ - SUBJECT_TYPE --to-format=json"
    command = command.split(" ")
    subp = subprocess.run(command, capture_output=True, encoding="UTF-8")
    assert subp.returncode == 0
    assert "{'SUBJECT_TYPE': {'Human'}}" == subp.stdout[:-1]


def test_extract_metadata_error_recovery(teardown_module):
    """Test that the extract metadata command can recover from an error."""
    command = "python -m mwtab extract metadata tests/example_data/files_to_test_error_recovery tests/example_data/tmp/ SUBJECT_TYPE"
    command = command.split(" ")
    subp = subprocess.run(command, capture_output=True, encoding="UTF-8")
    assert subp.returncode == 0
    assert "1_no_error.json" in subp.stdout
    assert "2_error.json" in subp.stdout
    assert "3_no_error.txt" in subp.stdout
    assert "ValueError: Blank input string retrieved from source." in subp.stdout


def test_extract_metadata_no_data_extracted(teardown_module):
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
def test_extract_metabolites_command(from_path, to_path, key, value, to_format, no_header, teardown_module):
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


def test_extract_metabolites_stdout():
    command = "python -m mwtab extract metabolites tests/example_data/mwtab_files/ - SU:SUBJECT_TYPE Human --to-format=json"
    command = command.split(" ")
    subp = subprocess.run(command, capture_output=True, encoding="UTF-8")
    assert subp.returncode == 0
    assert '{\n    "17-hydroxypregnenolone": {\n        "ST000122": {' in subp.stdout


def test_extract_metabolites_error_recovery(teardown_module):
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


def test_extract_metabolites_no_data_extracted(teardown_module):
    """Test that the extract metabolites command prints a message when no data is extracted and no file is made."""
    command = "python -m mwtab extract metabolites tests/example_data/mwtab_files/ tests/example_data/tmp/ SU:SUBJECT_TYPE Plant"
    command = command.split(" ")
    subp = subprocess.run(command, capture_output=True, encoding="UTF-8")
    assert not os.path.exists("tests/example_data/tmp/")
    assert subp.returncode == 0
    assert ("No metabolites extracted. No file was saved. " 
            "This is likely due to key value pairs filtering all of the studies out. "
            "Check your key value pairs and data.") in subp.stdout


