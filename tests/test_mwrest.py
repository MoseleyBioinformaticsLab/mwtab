import pytest
from mwtab.mwrest import GenericMWURL


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
])
def test_mwrest(kwds):
    test_mwurl = GenericMWURL(**kwds)
    assert test_mwurl.url == "https://www.metabolomicsworkbench.org/rest/{}/{}/{}/{}/{}".format(
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
])
def test_failure(kwds):
    expected_error = ValueError("Invalid Metabolomics Workbench analysis ID for a study (AN<6-digit integer>)")
    try:
        GenericMWURL(**kwds)
    except ValueError as e:
        assert e.args == expected_error.args
