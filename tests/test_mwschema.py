
from mwtab import mwschema


"""
Most functionality is tested through its utilization throughout the package, so these tests are just to get 100% coverage.
"""

def test_create_unit_error_message():
    assert ' should be a unitless integer. Ignore this when more complicated descriptions are required.' == mwschema.create_unit_error_message(integer=True, no_units=True)



