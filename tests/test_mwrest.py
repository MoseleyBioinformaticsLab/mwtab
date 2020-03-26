import pytest
from mwtab.mwrest import GenericMWURL, _pull_study_analysis


def test_study_analysis():
    test = _pull_study_analysis()
    assert test


@pytest.mark.parametrize("kwds", [
    ({"context": "study",
      "input_item": "analysis_id",
      "input_value": "AN000002",
      "output_item": "mwtab",
      'output_format': "txt"}),
    ({"context": "study",
      "input_item": "study_id",
      "input_value": "ST000001",
      "output_item": "mwtab",
      'output_format': "txt"}),
    ({"base_url": "https://www.test.org/rest/",
      "context": "study",
      "input_item": "study_id",
      "input_value": "ST000001",
      "output_item": "mwtab",
      'output_format': "txt"}),
])
def test_mwrest(kwds):
    test_mwurl = GenericMWURL(kwds)
    assert test_mwurl.url == test_mwurl.base_url + "/".join([
        kwds["context"],
        kwds["input_item"],
        kwds["input_value"],
        kwds["output_item"],
        kwds.get("output_format") or ""
    ])


@pytest.mark.parametrize("kwds", [
    ({"context": "study",
      "input_item": "analysis_id",
      "input_value": "ST000001",
      "output_item": "mwtab",
      'output_format': "txt"}),
    ({"context": "moverz",
      "input_item": "LIPIDS",
      "m/z_value": 49,
      "ion_type_value": "M+H",
      "m/z_tolerance_value": 0.1,
      'output_format': "txt"}),
    ({"context": "exactmass",
      "LIPID_abbreviation": "Test",
      "ion_type_value": "M+H"}),
])
def test_fail_mwrest(kwds):
    try:
        test_mwurl = GenericMWURL(kwds)
        assert False
    except Exception as e:
        assert type(e) == ValueError
