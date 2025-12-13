# -*- coding: utf-8 -*-


import pytest

from mwtab import tokenizer



def test_errors():
    with open("tests/example_data/other_mwtab_files/ST000022_AN000041_malformed_factor.txt", 'r') as f:
        text = f.read()
    
    with pytest.raises(ValueError, match = r".*Either a bar \('\| '\) separating 2 items.*") as exc_info:
        for token in tokenizer.tokenizer(text):
            pass
    
    err = exc_info.value
    assert hasattr(err, '__cause__')
    assert isinstance(err.__cause__, ValueError)
    assert err.__cause__.args[0] == ("Either a bar ('| ') separating 2 items is missing or "
                                     "there is an extra colon (':') in the factor "
                                     "key value pair, 'Disease Status::Cases'")
    
    
    with open("tests/example_data/other_mwtab_files/ST000022_AN000041_malformed_additional_data.txt", 'r') as f:
        text = f.read()
    
    with pytest.raises(ValueError, match = r".*Either a semicolon \('; '\) separating 2 items.*") as exc_info:
        for token in tokenizer.tokenizer(text):
            pass
    
    err = exc_info.value
    assert hasattr(err, '__cause__')
    assert isinstance(err.__cause__, ValueError)
    assert err.__cause__.args[0] == ("Either a semicolon ('; ') separating 2 items "
                                     "is missing or there is an extra equal sign "
                                     "('=') in the additional data key value pair, 'key1==value1'")
    
    
    with open("tests/example_data/other_mwtab_files/ST000022_AN000041_missing_tab.txt", 'r') as f:
        text = f.read()
    
    with pytest.raises(ValueError, match = r".*Expected a tab in the line.*") as exc_info:
        for token in tokenizer.tokenizer(text):
            pass
    
    err = exc_info.value
    assert hasattr(err, '__cause__')
    assert isinstance(err.__cause__, ValueError)
    assert err.__cause__.args[0] == "Expected a tab in the line."
    
    
    with open("tests/example_data/other_mwtab_files/ST000022_AN000041_malformed_SSF_line.txt", 'r') as f:
        text = f.read()
    
    with pytest.raises(IndexError, match = r".*LINE WITH ERROR:\n\t'SUBJECT_SAMPLE_FACTORS           \\t-C0559Disease Status:Cases\\t'") as exc_info:
        for token in tokenizer.tokenizer(text):
            pass
    




