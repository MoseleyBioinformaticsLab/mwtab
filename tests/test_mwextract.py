import pytest

from mwtab import mwextract
import mwtab


"""
Most functionality is tested through its utilization throughout the package, so these tests are just to get 100% coverage.
"""



def test_SetEncoder():
    with pytest.raises(TypeError, match = r'Object of type str is not JSON serializable'):
        mwextract.SetEncoder().default('a')

def test_extract_metabolites_exception():
    mwtabfile = mwtab.mwtab.MWTabFile('asdf') 
    mwtabfile['MS_METABOLITE_DATA'] = {'Data':[{'key1':'asdf'}]}
    assert mwextract.extract_metabolites([(mwtabfile, None)], [lambda x: True]) == {}
    









