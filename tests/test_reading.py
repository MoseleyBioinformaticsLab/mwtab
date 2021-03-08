import pytest
import mwtab


@pytest.mark.parametrize("files_source", [
    "204",
    "AN000204",
    "https://www.metabolomicsworkbench.org/rest/study/analysis_id/AN000204/mwtab/txt",
    "tests/example_data/mwtab_files/ST000122_AN000204.txt",
    "tests/example_data/mwtab_files/ST000122_AN000204.json",
])
def test_single_file_reading(files_source):
    mwtabfile_generator = mwtab.read_files(files_source)
    mwtabfile = next(mwtabfile_generator)
    assert mwtabfile.study_id == "ST000122"
    assert mwtabfile.analysis_id == "AN000204"


@pytest.mark.parametrize("files_source", [
    "tests/example_data/mwtab_files",
    "tests/example_data/mwtab_files.zip",
    "tests/example_data/mwtab_files.tar.gz",
    "tests/example_data/mwtab_files.tar.bz2"
])
def test_multiple_reading(files_source):
    mwtabfile_generator = mwtab.read_files(files_source)
    mwtabfiles_list = list(mwtabfile_generator)
    mwtabfiles_study_ids_set = set(mwf.study_id for mwf in mwtabfiles_list)
    mwtabfiles_analysis_ids_set = set(mwf.analysis_id for mwf in mwtabfiles_list)
    assert mwtabfiles_study_ids_set == {"ST000122"}
    assert mwtabfiles_analysis_ids_set == {"AN000204"}
