import os
import pytest
import mwtab
from json import loads

from fixtures import teardown_module

from mwtab.converter import Converter


teardown_module = teardown_module

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
        print(key)
        assert dict1[key] == dict2[key]


@pytest.mark.parametrize("mwtab_file_path, json_file_path", [
    ("tests/example_data/mwtab_files/ST000122_AN000204.txt", "tests/example_data/mwtab_files/ST000122_AN000204.json")
])
def test_convert_mwtab_to_json(mwtab_file_path, json_file_path, teardown_module):
    """

    """
    # convert given mwTab file to JSON
    mwfile = next(mwtab.read_files(mwtab_file_path))
    if not os.path.exists("tests/example_data/tmp/"):
        os.makedirs("tests/example_data/tmp/")
    with open("tests/example_data/tmp/tmp.json", "w", encoding="utf-8") as f:
        mwfile.write(f, file_format="json")
        f.close()

    # open files
    with open("tests/example_data/tmp/tmp.json", "r", encoding="utf-8") as f:
        mwtab_file = loads(f.read())
    with open(json_file_path, "r", encoding="utf-8") as f:
        json_file = loads(f.read())

    # assert both files contain the same sections
    assert not set(mwtab_file.keys()) ^ set(json_file.keys())

    # Assert item sections are equal
    for section_key in ITEM_SECTIONS:
        print(section_key)
        if section_key in set(mwtab_file.keys()) & set(json_file.keys()):
            compare_item_sections(mwtab_file[section_key], json_file[section_key])

    # assert MS_METABOLITE_DATA or NMR_METABOLITE_DATA sections are the same


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
    assert mwtabfiles_study_ids_set.issubset({"ST000122"})
    assert mwtabfiles_analysis_ids_set.issubset({"AN000204"})


@pytest.mark.parametrize("from_path, to_path, from_format, to_format, error, message", [
    ("nonexistent_file.txt", "tests/example_data/tmp/json/ST000122_AN000204.json", "mwtab", "json", FileNotFoundError, r'No such file or directory: "nonexistent_file.txt"'),
    ("tests/example_data/other_test_data/unknown_file_format.asdf", "tests/example_data/tmp/json/ST000122_AN000204.json", "mwtab", "json", TypeError, r'Unknown input file format: "tests/example_data/other_test_data/unknown_file_format.asdf"'),
])
def test_convert_exceptions(teardown_module, from_path, to_path, from_format, to_format, error, message, mocker):
    converter = Converter(from_path=from_path,
                          to_path=to_path,
                          from_format=from_format,
                          to_format=to_format)
    
    # I can't find a way to trigger the last else without mocking, but I don't want to remove the line either.
    if error is TypeError:
        mocker.patch('mwtab.converter.os.path.isfile', side_effect = [False])
    
    with pytest.raises(error, match = message):
        converter.convert()


@pytest.mark.parametrize("from_path, to_path, from_format, to_format, error, message", [
    ("tests/example_data/mwtab_files/ST000122_AN000204.txt", "bad_path.gz", "mwtab", "json", TypeError, r'Many-to-one conversion, cannot convert "tests/example_data/mwtab_files/ST000122_AN000204.txt" into "bad_path.gz"'),
    ("tests/example_data/mwtab_files/ST000122_AN000204.txt", "tests/example_data/other_test_data/unknown_file_format.asdf", "mwtab", "json", TypeError, r'Unknown output file format: "tests/example_data/other_test_data/unknown_file_format.asdf"'),
])
def test_many_to_many_exceptions(teardown_module, from_path, to_path, from_format, to_format, error, message):
    converter = Converter(from_path=from_path,
                          to_path=to_path,
                          from_format=from_format,
                          to_format=to_format)
    
    # Can't get the last else line to trigger any other way than to do something like this.
    if 'unknown_file_format' in to_path:
        converter.file_generator.to_path_compression = '.asdf'
    
    with pytest.raises(error, match = message):
        converter._many_to_many()


@pytest.mark.parametrize("from_path, to_path, from_format, to_format, error, message", [
    ("tests/example_data/mwtab_files/ST000122_AN000204.txt", "bad_path.tar", "mwtab", "json", TypeError, r'One-to-many conversion, cannot convert "tests/example_data/mwtab_files/ST000122_AN000204.txt" into "bad_path.tar"'),
    ("tests/example_data/mwtab_files/ST000122_AN000204.txt", "tests/example_data/other_test_data/unknown_file_format.asdf", "mwtab", "json", TypeError, r'Unknown output file format: "tests/example_data/other_test_data/unknown_file_format.asdf"'),
])
def test_one_to_one_exceptions(teardown_module, from_path, to_path, from_format, to_format, error, message):
    converter = Converter(from_path=from_path,
                          to_path=to_path,
                          from_format=from_format,
                          to_format=to_format)
    
    # Can't get the last else line to trigger any other way than to do something like this.
    if 'unknown_file_format' in to_path:
        converter.file_generator.to_path_compression = '.asdf'
    
    with pytest.raises(error, match = message):
        converter._one_to_one()


def test_to_dir_exceptions(teardown_module, capsys, mocker):
    converter = Converter("tests/example_data/mwtab_files/ST000122_AN000204.txt", "tests/example_data/tmp/json/ST000122_AN000204.json", "mwtab", "json")
    
    # Simulate an exception with this mock.
    mocker.patch('mwtab.converter.os.path.exists', side_effect = [True, False])
    converter._to_dir(converter.file_generator)
    captured = capsys.readouterr()
    assert 'Something went wrong when trying to convert' in captured.out


def test_to_zipfile_exceptions(teardown_module, capsys, mocker):
    converter = Converter("tests/example_data/mwtab_files/ST000122_AN000204.txt", "tests/example_data/tmp/json.zip", "mwtab", "json")
    
    # Simulate an exception with this mock.
    mocker.patch('mwtab.converter.Converter._output_path', side_effect = [TypeError])
    
    # None of the private methods check for this, since convert() does it, so it must be created here.
    if not os.path.exists('tests/example_data/tmp'):
        os.makedirs('tests/example_data/tmp')
    
    converter._to_zipfile(converter.file_generator)
    captured = capsys.readouterr()
    assert 'Something went wrong when trying to convert' in captured.out


def test_to_tarfile_exceptions(teardown_module, capsys, mocker):
    converter = Converter("tests/example_data/mwtab_files/ST000122_AN000204.txt", "tests/example_data/tmp/json.tar", "mwtab", "json")
    
    # This is not necessary for the exception, but tests the last else line.
    converter.file_generator.to_path_compression = 'asdf'
    
    # Simulate an exception with this mock.
    mocker.patch('mwtab.converter.Converter._output_path', side_effect = [TypeError])
    
    # None of the private methods check for this, since convert() does it, so it must be created here.
    if not os.path.exists('tests/example_data/tmp'):
        os.makedirs('tests/example_data/tmp')
    
    converter._to_tarfile(converter.file_generator)
    captured = capsys.readouterr()
    assert 'Something went wrong when trying to convert' in captured.out


def test_to_bz2file_exceptions(teardown_module, capsys):
    converter = Converter("tests/example_data/mwtab_files/ST000122_AN000204.txt", "tests/example_data/tmp/json.bz2", "mwtab", "json")
    
    converter.file_generator.to_format = 0
    
    # None of the private methods check for this, since convert() does it, so it must be created here.
    if not os.path.exists('tests/example_data/tmp'):
        os.makedirs('tests/example_data/tmp')
    
    converter._to_bz2file(converter.file_generator)
    captured = capsys.readouterr()
    assert 'Something went wrong when trying to convert' in captured.out


def test_to_gzipfile_exceptions(teardown_module, capsys):
    converter = Converter("tests/example_data/mwtab_files/ST000122_AN000204.txt", "tests/example_data/tmp/json.gz", "mwtab", "json")
    
    converter.file_generator.to_format = 0
    
    # None of the private methods check for this, since convert() does it, so it must be created here.
    if not os.path.exists('tests/example_data/tmp'):
        os.makedirs('tests/example_data/tmp')
    
    converter._to_gzipfile(converter.file_generator)
    captured = capsys.readouterr()
    assert 'Something went wrong when trying to convert' in captured.out


def test_to_textfile_exceptions(teardown_module, capsys, mocker):
    converter = Converter("tests/example_data/mwtab_files/ST000122_AN000204.txt", "tests/example_data/tmp/txt.txt", "mwtab", "mwtab")
    
    mocker.patch('mwtab.fileio.mwtab.MWTabFile.writestr', side_effect = [TypeError])
    
    # None of the private methods check for this, since convert() does it, so it must be created here.
    if not os.path.exists('tests/example_data/tmp'):
        os.makedirs('tests/example_data/tmp')
    
    converter._to_textfile(converter.file_generator)
    captured = capsys.readouterr()
    assert 'Something went wrong when trying to convert' in captured.out




