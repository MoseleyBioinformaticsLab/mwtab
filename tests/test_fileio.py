# -*- coding: utf-8 -*-

import pytest

from mwtab import fileio


"""
Most functionality is tested through its utilization throughout the package, so these tests are just to get 100% coverage.
"""

def test_return_correct_yield():
    with pytest.raises(TypeError):
        fileio._return_correct_yield({}, TypeError(), False)


def test_generate_filenames_exception(mocker):
    mocker.patch('mwtab.fileio.os.path.isdir', side_effect = [Exception])
    with pytest.raises(Exception):
        next(fileio._generate_filenames(['tests/example_data/mwtab_files/']))


def test_generate_handles_exception():
    generator = fileio._generate_handles([('name1', Exception()), ('name2', None)])
    with pytest.raises(Exception):
        next(generator)
    with pytest.raises(StopIteration):
        next(generator)

def test_generate_handles_exception2():
    generator = fileio._generate_handles([('name1', Exception()), ('name2', None)], True)
    value1 = next(generator)
    assert value1[0] is None
    assert value1[1] == 'name1'
    assert isinstance(value1[2], Exception)
    
    value2 = next(generator)
    assert value2[0] is None
    assert value2[1] == 'name2'
    assert isinstance(value2[2], FileNotFoundError)


def test_read_files_exceptions(mocker):
    with pytest.raises(TypeError, match = r'Unknown file source.'):
        next(fileio.read_files(['some_path']))
    
    mocker.patch('mwtab.fileio._generate_filenames', side_effect = [Exception])
    with pytest.raises(Exception):
        next(fileio.read_files(['some_path']))

def test_read_files_exceptions2(mocker, capsys):
    fileio.VERBOSE = True
    mocker.patch('mwtab.mwtab.MWTabFile.read', side_effect = [Exception])
    with pytest.raises(Exception):
        next(fileio.read_files(['tests/example_data/mwtab_files/ST000122_AN000204.txt']))
    captured = capsys.readouterr()
    assert 'Error processing file:' in captured.out


def test_GenericFilePath():
    web_address = 'https://github.com/MoseleyBioinformaticsLab/mwtab/archive/refs/heads/main.zip'
    filehandle, source = next(fileio.GenericFilePath(web_address).open())
    filehandle.read()
    assert web_address in source

def test_GenericFilePath_is_url(mocker):
    mocker.patch('mwtab.fileio.urlparse', side_effect = [ValueError])
    assert fileio.GenericFilePath.is_url('asdf') == False


def test_read_lines(mocker, capsys):
    """We are really testing the ReadLines class, but doing it through the read_lines function."""
    fileio.VERBOSE = True
    lines = next(fileio.read_lines('tests/example_data/mwtab_files/ST000122_AN000204.txt')).lines
    captured = capsys.readouterr()
    assert 'Processed file: ' in captured.out
    assert lines[0].startswith('#METABOLOMICS WORKBENCH')
        
    mocker.patch('mwtab.fileio._generate_handles', side_effect = [((open('tests/example_data/mwtab_files/ST000122_AN000204.txt', 'rb'), 'asdf', None) for a in [1])])
    lines = next(fileio.read_lines('tests/example_data/mwtab_files/ST000122_AN000204.txt')).lines
    assert lines[0].startswith('#METABOLOMICS WORKBENCH')

def test_read_lines_exceptions(mocker, capsys):
    """We are really testing the ReadLines class, but doing it through the read_lines function."""
    fileio.VERBOSE = True
    class mock_class:
        def __init__(self, *args, **kwargs):
            pass
        def read(self):
            return 1
        def close(self):
            pass
    mocker.patch('mwtab.fileio._generate_handles', side_effect = [((mock_class(), 'asdf', None) for a in [1])])
    with pytest.raises(TypeError, match = r"Expecting <class 'str'> or <class 'bytes'>,"):
        next(fileio.read_lines(['tests/example_data/mwtab_files/ST000122_AN000204.txt']))
    captured = capsys.readouterr()
    assert 'Error processing file:' in captured.out

