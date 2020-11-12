import os
import shutil
import pytest
import mwtab
from json import loads
from mwtab.converter import Converter


ITEM_SECTIONS = {
    # "METABOLOMICS WORKBENCH",
    "PROJECT",
    "STUDY",
    "ANALYSIS",
    "SUBJECT",
    "COLLECTION",
    "TREATMENT",
    "SAMPLEPREP",
    "CHROMATOGRAPHY",
    "MS",
    "NMR",
}


def teardown_module(module):
    if os.path.exists("tests/example_data/tmp"):
        shutil.rmtree("tests/example_data/tmp")


def compare_item_sections(dict1, dict2):
    """
    Method for comparing the item sections of two given dictionaries.

    Helper method which asserts two item sections (dictionaries), section which only contain key-value item pairs, from
    two different `~mwtab.mwtab.MWTabFile` objects are equal.

    :param dict1: First dictionary representing mwTab file section containing key-value item pairs.
    :type dict1: :py:class:`collections.OrderedDict` or :py:class:`dict`
    :param dict2: Second dictionary representing mwTab file section containing key-value item pairs.
    :type dict2: :py:class:`collections.OrderedDict` or :py:class:`dict`
    """
    keys1 = set(dict1.keys())
    keys2 = set(dict2.keys())

    assert not keys1 ^ keys2

    for key in keys1 & keys2:
        assert dict1[key] == dict2[key]


@pytest.mark.parametrize("mwtab_file_path, json_file_path", [
    ("tests/example_data/converter_files/ST000001_AN000001.txt", "tests/example_data/converter_files/ST000001_AN000001.json"),
    ("tests/example_data/converter_files/ST001468_AN002444.txt", "tests/example_data/converter_files/ST001468_AN002444.json")
])
def test_convert_mwtab_to_json(mwtab_file_path, json_file_path):
    """

    """
    # convert given mwTab file to JSON
    mwfile = next(mwtab.read_files(mwtab_file_path))
    if not os.path.exists("tests/example_data/tmp/"):
        os.makedirs("tests/example_data/tmp/")
    with open("tests/example_data/tmp/tmp.json", "w") as f:
        mwfile.write(f, file_format="json")
        f.close()

    # open files
    with open("tests/example_data/tmp/tmp.json", "r") as f:
        mwtab_file = loads(f.read())
    with open(json_file_path, "r") as f:
        json_file = loads(f.read())

    # assert both files contain the same sections
    assert not set(mwtab_file.keys()) ^ set(json_file.keys())

    # Assert item sections are equal
    for section_key in ITEM_SECTIONS:
        if section_key in set(mwtab_file.keys()) & set(json_file.keys()):
            compare_item_sections(mwtab_file[section_key], json_file[section_key])

    # assert MS_METABOLITE_DATA or NMR_METABOLITE_DATA sections are the same


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
def test_converter_module(from_path, to_path, from_format, to_format):
    converter = Converter(from_path=from_path,
                          to_path=to_path,
                          from_format=from_format,
                          to_format=to_format)
    converter.convert()

    mwtabfile_generator = mwtab.read_files(to_path)
    mwtabfiles_list = list(mwtabfile_generator)
    mwtabfiles_study_ids_set = set(mwf.study_id for mwf in mwtabfiles_list)
    mwtabfiles_analysis_ids_set = set(mwf.analysis_id for mwf in mwtabfiles_list)
    assert mwtabfiles_study_ids_set.issubset({"ST000001", "ST000002"})
    assert mwtabfiles_analysis_ids_set.issubset({"AN000001", "AN000002"})


if __name__ == "__main__":
    test_convert_mwtab_to_json(
        "tests/example_data/converter_files/ST000001_AN000001.txt",
        "tests/example_data/converter_files/ST000001_AN000001.json"
    )
