import pytest
from mwtab.mwrest import GenericMWURL, _pull_study_analysis


def test_study_analysis():
    test = _pull_study_analysis()
    assert test


@pytest.mark.parametrize("kwds", [
    ({'context': 'study',
      'input item': 'analysis_id',
      'input value': "AN000002",
      'output item': 'mwtab',
      'output format': "txt"}),
    ({'context': 'study',
      'input item': 'study_id',
      'input value': "ST000001",
      'output item': 'mwtab'}),
    ({'base url': "https://www.test.org/rest/",
      'context': 'study',
      'input item': 'study_id',
      'input value': "ST000001",
      'output item': 'mwtab'})
])
def test_mwrest(kwds):
    test_mwurl = GenericMWURL(**kwds)
    assert test_mwurl.url == test_mwurl.base_url + "{}/{}/{}/{}/{}".format(
        kwds["context"],
        kwds["input item"],
        kwds["input value"],
        kwds["output item"],
        kwds.get("output format") or ""
    )


@pytest.mark.parametrize("kwds", [
    ({'context': 'study',
      'input item': 'analysis_id',
      'input value': "ST000002",
      'output item': 'mwtab',
      'output format': "txt"}),
    ({'context': 'moverz',
      'input item': 'LIPIDS',
      'm/z value': 49,
      'ion type value': "M+H",
      'm/z tolerance value': 0.1,
      'output format': "txt"}),
    ({'context': 'exactmass',
      'LIPID abbreviation': 'Test',
      'ion type value': "M+H"}),
])
def test_fail_mwrest(kwds):
    try:
        test_mwurl = GenericMWURL(**kwds)
        assert False
    except Exception as e:
        assert type(e) == ValueError
