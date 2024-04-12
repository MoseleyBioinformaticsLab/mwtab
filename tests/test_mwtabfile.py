# -*- coding: utf-8 -*-

import mwtab
import os
import shutil
from json import loads


def teardown_module(module):
    if os.path.exists("tests/example_data/tmp"):
        shutil.rmtree("tests/example_data/tmp")
        

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

