# -*- coding: utf-8 -*-

import pandas

from mwtab import metadata_column_matching


"""
Most functionality is tested through its utilization throughout the package, so these tests are just to get 100% coverage.
"""


def test_ValueMatcher():
    vm1 = metadata_column_matching.ValueMatcher('integer', None, None)
    series = pandas.Series(['a', '1', '2.3', '4'])
    series_match = vm1.series_match(series, None)
    expected_match = pandas.Series([False, True, False, True])
    assert series_match.equals(expected_match)
    
    # These conflict, so nothing should match.
    vm1 = metadata_column_matching.ValueMatcher('integer', r'a', None)
    series = pandas.Series(['a', '1', '2.3', '4'])
    series_match = vm1.series_match(series, None)
    expected_match = pandas.Series([False, False, False, False])
    assert series_match.equals(expected_match)
    
    vm1 = metadata_column_matching.ValueMatcher('numeric', None, None)
    series = pandas.Series(['a', '1', '2.3', '4'])
    series_match = vm1.series_match(series, None)
    expected_match = pandas.Series([False, True, True, True])
    assert series_match.equals(expected_match)












