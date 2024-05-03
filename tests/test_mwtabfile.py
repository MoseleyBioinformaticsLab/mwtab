# -*- coding: utf-8 -*-

import mwtab
import os
import shutil
from json import loads
import time
import pathlib


def teardown_module(module):
    path = pathlib.Path("tests/example_data/tmp/")
    if os.path.exists(path):
        shutil.rmtree(path)
        time_to_wait=10
        time_counter = 0
        while path.exists():
            time.sleep(1)
            time_counter += 1
            if time_counter > time_to_wait:
                raise FileExistsError(path + " was not deleted within " + str(time_to_wait) + " seconds, so it is assumed that it won't be and something went wrong.")

        

def test_keys_reorder():
    """Test that the keys are reoodered to match what the Workbench wants."""
    test_dict = {"MS_METABOLITE_DATA": {"Data":[{"attribute1":"qwer", "attribute2":"zxcv", "Metabolite":"asdf"}], "Units":"some unit"}, 
                  "ANALYSIS":{"asdf":"asdf"},
                  "SUBJECT_SAMPLE_FACTORS":[
                                        {
                                        "Factors":{"Feeeding":"Ad lib","Running Capacity":"High"},
                                        "Subject ID":"-",
                                        "Sample ID":"S00009477"
                                        },
                                        {
                                        "Subject ID":"-",
                                        "Sample ID":"S00009478",
                                        "Factors":{"Feeeding":"Ad lib","Running Capacity":"High"}
                                        }]}
    
    expected_dict = {'ANALYSIS': {'asdf': 'asdf'},
                      'SUBJECT_SAMPLE_FACTORS': [
                          {'Subject ID': '-',
                          'Sample ID': 'S00009477',
                          'Factors':
                                {'Feeeding': 'Ad lib',
                                'Running Capacity': 'High'}},
                          {'Subject ID': '-',
                          'Sample ID': 'S00009478',
                          'Factors':
                                {'Feeeding': 'Ad lib',
                                'Running Capacity': 'High'}}],
                      'MS_METABOLITE_DATA':{
                                          'Units': 'some unit',
                                          'Data':[
                                                  {'Metabolite': 'asdf',
                                                    'attribute1': 'qwer',
                                                    'attribute2': 'zxcv'}]}}

    mwtabfile = mwtab.mwtab.MWTabFile("some file path")

    for key, value in test_dict.items():
        mwtabfile[key] = value

    if not os.path.exists("tests/example_data/tmp/"):
        os.makedirs("tests/example_data/tmp/")
        
    with open("tests/example_data/tmp/tmp.json", "w", encoding="utf-8") as f:
        mwtabfile.write(f, file_format="json")
    
    with open("tests/example_data/tmp/tmp.json", "r", encoding="utf-8") as f:
        json_file = loads(f.read())
        
    assert expected_dict == json_file


def test_read_in_and_reorder_keys():
    """Test that a file with incorrect key order gets reordered."""
    
    if not os.path.exists("tests/example_data/tmp/"):
        os.makedirs("tests/example_data/tmp/")
    
    mwtabfile = mwtab.mwtab.MWTabFile("tests/example_data/other_mwtab_files/incorrect_section_order.json")
        
    with open("tests/example_data/other_mwtab_files/incorrect_section_order.json", "r", encoding="utf-8") as f:
        mwtabfile.read(f)
    
    with open("tests/example_data/tmp/tmp.json", "w", encoding="utf-8") as f:
        mwtabfile.write(f, file_format="json")
    
    with open("tests/example_data/tmp/tmp.json", "r", encoding="utf-8") as f:
        json_file = loads(f.read())
    
    with open("tests/example_data/other_mwtab_files/corrected_section_order.json", "r", encoding="utf-8") as f:
        expected_dict = loads(f.read())
        
    assert expected_dict == json_file


def test_read_in_duplicate_keys_json():
    """Test that a file with duplicate keys is handled correctly for JSON files."""
    
    if not os.path.exists("tests/example_data/tmp/"):
        os.makedirs("tests/example_data/tmp/")
    
    mwtabfile = mwtab.mwtab.MWTabFile("tests/example_data/other_mwtab_files/ST000122_AN000204_duplicate_keys.json")
        
    with open("tests/example_data/other_mwtab_files/ST000122_AN000204_duplicate_keys.json", "r", encoding="utf-8") as f:
        mwtabfile.read(f)
    
    assert isinstance(mwtabfile["SUBJECT_SAMPLE_FACTORS"][0]["Additional sample data"]['key_1'], mwtab.mwtab._duplicate_key_list)
    
    with open("tests/example_data/tmp/tmp.json", "w", encoding="utf-8") as f:
        mwtabfile.write(f, file_format="json")
    
    new_mwtabfile = mwtab.mwtab.MWTabFile("tests/example_data/tmp/tmp.json")
    with open("tests/example_data/tmp/tmp.json", "r", encoding="utf-8") as f:
        new_mwtabfile.read(f)
    
    assert isinstance(new_mwtabfile["SUBJECT_SAMPLE_FACTORS"][0]["Additional sample data"]['key_1'], mwtab.mwtab._duplicate_key_list)


def test_read_in_duplicate_keys_tab():
    """Test that a file with duplicate keys is handled correctly for tab files."""
    
    if not os.path.exists("tests/example_data/tmp/"):
        os.makedirs("tests/example_data/tmp/")
    
    mwtabfile = mwtab.mwtab.MWTabFile("tests/example_data/other_mwtab_files/ST000122_AN000204_duplicate_keys.txt")
        
    with open("tests/example_data/other_mwtab_files/ST000122_AN000204_duplicate_keys.txt", "r", encoding="utf-8") as f:
        mwtabfile.read(f)
    
    assert isinstance(mwtabfile["SUBJECT_SAMPLE_FACTORS"][0]["Additional sample data"]['key_1'], mwtab.mwtab._duplicate_key_list)
    
    with open("tests/example_data/tmp/tmp.txt", "w", encoding="utf-8") as f:
        mwtabfile.write(f, file_format="mwtab")
    
    new_mwtabfile = mwtab.mwtab.MWTabFile("tests/example_data/tmp/tmp.txt")
    with open("tests/example_data/tmp/tmp.txt", "r", encoding="utf-8") as f:
        new_mwtabfile.read(f)
    
    assert isinstance(new_mwtabfile["SUBJECT_SAMPLE_FACTORS"][0]["Additional sample data"]['key_1'], mwtab.mwtab._duplicate_key_list)    


def test_validate():
    """Test that the validate method validates the object."""
    
    mwtabfile = mwtab.mwtab.MWTabFile("tests/example_data/other_mwtab_files/ST000122_AN000204_duplicate_keys.txt")
        
    with open("tests/example_data/other_mwtab_files/ST000122_AN000204_duplicate_keys.txt", "r", encoding="utf-8") as f:
        mwtabfile.read(f)
    
    _, errors = mwtabfile.validate(verbose=False)
    
    assert "duplicate keys" in errors
    
    
def test_from_dict():
    """Test that the from_dict method works to create a new MWTabFile object."""
    
    with open("tests/example_data/other_mwtab_files/incorrect_section_order.json", "r", encoding="utf-8") as f:
        json_file = loads(f.read())
    
    mwtabfile = mwtab.mwtab.MWTabFile.from_dict(json_file)
    
    assert mwtabfile.study_id == "ST000000"


def test_properties():
    """Test that the study_id, analysis_id, and header properties behave as expected."""
    
    mwtabfile = mwtab.mwtab.MWTabFile("tests/example_data/other_mwtab_files/ST000122_AN000204_duplicate_keys.txt")
        
    with open("tests/example_data/other_mwtab_files/ST000122_AN000204_duplicate_keys.txt", "r", encoding="utf-8") as f:
        mwtabfile.read(f)
    
    assert mwtabfile.study_id == "ST000122"
    assert mwtabfile.analysis_id == "AN000204"
    assert mwtabfile.header == "#METABOLOMICS WORKBENCH STUDY_ID:ST000122 ANALYSIS_ID:AN000204 PROJECT_ID:PR000109"
    
    temp = mwtabfile["METABOLOMICS WORKBENCH"]
    del mwtabfile["METABOLOMICS WORKBENCH"]
    
    assert mwtabfile.study_id is None
    assert mwtabfile.analysis_id is None
    assert mwtabfile.header is None
    
    mwtabfile["METABOLOMICS WORKBENCH"] = temp
    
    assert mwtabfile.study_id == "ST000122"
    assert mwtabfile.analysis_id == "AN000204"
    assert mwtabfile.header == "#METABOLOMICS WORKBENCH STUDY_ID:ST000122 ANALYSIS_ID:AN000204 PROJECT_ID:PR000109"
    
    mwtabfile.study_id = "asdf"
    mwtabfile.analysis_id = "qwer"
    mwtabfile.header = "zxcv"
    
    assert mwtabfile.study_id == "asdf"
    assert mwtabfile.analysis_id == "qwer"
    assert mwtabfile.header == "zxcv"




