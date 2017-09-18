import pytest
import mwtab


def test_from_local_file():
    mwtabfile_generator = mwtab.read_files("tests/example_data/mwtab_files/ST000001_AN000001.txt",
                                           "tests/example_data/mwtab_files/ST000002_AN000002.txt")
    mwtabfile1 = next(mwtabfile_generator)
    mwtabfile2 = next(mwtabfile_generator)
    assert mwtabfile1.study_id == "ST000001" and mwtabfile2.study_id == "ST000002"
    assert mwtabfile1.analysis_id == "AN000001" and mwtabfile2.analysis_id == "AN000002"


def test_from_analysis_id():
    mwtabfile_generator = mwtab.read_files("1")
    mwtabfile = next(mwtabfile_generator)
    assert mwtabfile.study_id == "ST000001"
    assert mwtabfile.analysis_id == "AN000001"


@pytest.mark.parametrize("files_source", [
    "tests/example_data/mwtab_files",
    "tests/example_data/mwtab_files.zip",
    "tests/example_data/mwtab_files.tar.gz",
    "tests/example_data/mwtab_files.tar.bz2"
])
def test_reading(files_source):
    mwtabfile_generator = mwtab.read_files(files_source)
    mwtabfiles_list = list(mwtabfile_generator)
    mwtabfiles_study_ids_set = set(mwf.study_id for mwf in mwtabfiles_list)
    mwtabfiles_analysis_ids_set = set(mwf.analysis_id for mwf in mwtabfiles_list)
    assert mwtabfiles_study_ids_set == {"ST000001", "ST000002"}
    assert mwtabfiles_analysis_ids_set == {"AN000001", "AN000002"}
