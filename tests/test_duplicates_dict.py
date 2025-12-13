# -*- coding: utf-8 -*-

import copy

from mwtab.duplicates_dict import DuplicatesDict


"""
Most functionality is tested through its utilization throughout the package, so these tests are just to get 100% coverage.
"""

def test_init_from_dd():
    test_dict = DuplicatesDict()
    test_dict['a'] = 1
    test_dict2 = DuplicatesDict(test_dict)
    assert test_dict == test_dict2

def test_len():
    test_dict = DuplicatesDict()
    test_dict['a'] = 1
    assert len(test_dict) == 1

def test_delitem():
    test_dict = DuplicatesDict()
    test_dict['a'] = 1
    del test_dict['a']
    assert 'a' not in test_dict

def test_raw_keys():
    test_dict = DuplicatesDict()
    test_dict['a'] = 1
    test_dict['a'] = 2
    raw_keys = test_dict.raw_keys()
    assert 'a' in raw_keys
    assert 'a{{{_1_}}}' in raw_keys
    
def test_update():
    test_dict = DuplicatesDict()
    test_dict.update({'a':1})
    assert 'a' in test_dict
    assert test_dict['a'] == 1

def test_copy():
    test_dict = DuplicatesDict()
    test_dict['a'] = 1
    test_dict2 = copy.copy(test_dict)
    assert test_dict == test_dict2
    assert id(test_dict) != id(test_dict2)

def test_not_equal():
    test_dict = DuplicatesDict()
    test_dict['a'] = 1
    test_dict2 = DuplicatesDict()
    test_dict2['b'] = 2
    assert test_dict != test_dict2




