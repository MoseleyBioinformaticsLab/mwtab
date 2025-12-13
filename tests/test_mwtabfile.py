# -*- coding: utf-8 -*-

import mwtab
import os
from json import loads
import copy
import io
import json

import pandas
import pytest
# This is an autouse module, so it is being used by every test simply by importing without calling it directly.
from fixtures import teardown_module_auto, init_tmp_dir

init_tmp_dir = init_tmp_dir



def test_read():
    # mwtab
    mwtabfile = mwtab.mwtab.MWTabFile("tests/example_data/other_mwtab_files/ST000122_AN000204_duplicate_keys.txt")
    with open("tests/example_data/other_mwtab_files/ST000122_AN000204_duplicate_keys.txt", "r", encoding="utf-8") as f:
        mwtabfile.read(f)
    
    assert mwtabfile['MS_METABOLITE_DATA']['Data'][0]['Metabolite'] == '17-hydroxypregnenolone'
    assert mwtabfile._input_format == 'mwtab'
    
    # json
    mwtabfile2 = mwtab.mwtab.MWTabFile("tests/example_data/other_mwtab_files/ST000122_AN000204_duplicate_keys.json")
    with open("tests/example_data/other_mwtab_files/ST000122_AN000204_duplicate_keys.json", "r", encoding="utf-8") as f:
        mwtabfile2.read(f)
    
    assert mwtabfile2['MS_METABOLITE_DATA']['Data'][0]['Metabolite'] == '17-hydroxypregnenolone'
    assert mwtabfile2._input_format == 'json'
    
    # bytes
    mwtabfile2 = mwtab.mwtab.MWTabFile("tests/example_data/other_mwtab_files/ST000122_AN000204_duplicate_keys.json")
    with open("tests/example_data/other_mwtab_files/ST000122_AN000204_duplicate_keys.json", "rb") as f:
        mwtabfile2.read(f)
    
    assert mwtabfile2['MS_METABOLITE_DATA']['Data'][0]['Metabolite'] == '17-hydroxypregnenolone'
    assert mwtabfile2._input_format == 'json'
    
    # binned json
    mwtabfile = mwtab.mwtab.MWTabFile("tests/example_data/other_mwtab_files/ST000122_AN000204_binned.json")
    with open("tests/example_data/other_mwtab_files/ST000122_AN000204_binned.json", "r", encoding="utf-8") as f:
        mwtabfile.read(f)
    
    assert mwtabfile['NMR_BINNED_DATA']['Data'][0]['Metabolite'] == '17-hydroxypregnenolone'
    
    # binned mwtab
    mwtabfile = mwtab.mwtab.MWTabFile("tests/example_data/other_mwtab_files/ST000022_AN000041.txt")
    with open("tests/example_data/other_mwtab_files/ST000022_AN000041.txt", "r", encoding="utf-8") as f:
        mwtabfile.read(f)
    
    assert mwtabfile['NMR_BINNED_DATA']['Data'][0]['Metabolite'] == '0.4...0.46'
    
def test_read_errors():
    with pytest.raises(TypeError, match = r'^Unknown file format'):
        mwtabfile = mwtab.mwtab.MWTabFile("tests/example_data/other_mwtab_files/bad_file.txt")
        with open("tests/example_data/other_mwtab_files/bad_file.txt", "r", encoding="utf-8") as f:
            mwtabfile.read(f)
    
    with pytest.raises(KeyError, match = r'^\'Given JSON contains duplicate keys at the highest level\.'):
        mwtabfile = mwtab.mwtab.MWTabFile("tests/example_data/other_mwtab_files/ST000122_AN000204_duplicate_keys_highest_level.json", duplicate_keys=True)
        with open("tests/example_data/other_mwtab_files/ST000122_AN000204_duplicate_keys_highest_level.json", "r", encoding="utf-8") as f:
            mwtabfile.read(f)
    
    with pytest.raises(TypeError, match = r'^Given JSON contains non-dictionary values in \["MS_METABOLITE_DATA"\]\["Data"\]\.'):
        mwtabfile = mwtab.mwtab.MWTabFile("tests/example_data/other_mwtab_files/ST000122_AN000204_bad_table_type.json")
        with open("tests/example_data/other_mwtab_files/ST000122_AN000204_bad_table_type.json", "r", encoding="utf-8") as f:
            mwtabfile.read(f)

def test_read_with_force():
    mwtabfile = mwtab.mwtab.MWTabFile("tests/example_data/other_mwtab_files/ST000122_AN000204_bad_table_type.json", force=True)
    with open("tests/example_data/other_mwtab_files/ST000122_AN000204_bad_table_type.json", "r", encoding="utf-8") as f:
        mwtabfile.read(f)
    assert mwtabfile['MS_METABOLITE_DATA']['Data'][0] != 'asdf'
    

def test_reading_results_file():
    results_file_dict = {'filename': 'ST000071_AN000111_Results.txt',
                          'UNITS': 'Peak area',
                          'Has m/z': 'Yes',
                          'Has RT': 'Yes',
                          'RT units': 'Minutes'}
    mwtabfile = mwtab.mwtab.MWTabFile("tests/example_data/other_mwtab_files/ST000022_AN000041_results_file.txt")
    with open("tests/example_data/other_mwtab_files/ST000022_AN000041_results_file.txt", "r", encoding="utf-8") as f:
        mwtabfile.read(f)
    
    assert mwtabfile['NM']['NMR_RESULTS_FILE'] == results_file_dict
    
    
    mwtabfile = mwtab.mwtab.MWTabFile("tests/example_data/other_mwtab_files/ST000122_AN000204_results_file.txt")
    with open("tests/example_data/other_mwtab_files/ST000122_AN000204_results_file.txt", "r", encoding="utf-8") as f:
        mwtabfile.read(f)
    
    assert mwtabfile['MS']['MS_RESULTS_FILE'] == results_file_dict
    
    
    mwtabfile = mwtab.mwtab.MWTabFile("tests/example_data/other_mwtab_files/ST000022_AN000041_results_file_bad_location.txt")
    with open("tests/example_data/other_mwtab_files/ST000022_AN000041_results_file_bad_location.txt", "r", encoding="utf-8") as f:
        mwtabfile.read(f)
    
    assert mwtabfile['NM']['NMR_RESULTS_FILE'] == results_file_dict
    
    
    mwtabfile = mwtab.mwtab.MWTabFile("tests/example_data/other_mwtab_files/ST000122_AN000204_results_file_bad_location.txt")
    with open("tests/example_data/other_mwtab_files/ST000122_AN000204_results_file_bad_location.txt", "r", encoding="utf-8") as f:
        mwtabfile.read(f)
    
    assert mwtabfile['MS']['MS_RESULTS_FILE'] == results_file_dict
    

def test_reading_duplicate_subsections():
    mwtabfile = mwtab.mwtab.MWTabFile("tests/example_data/other_mwtab_files/ST000122_AN000204_duplicate_subsections.txt")
    with open("tests/example_data/other_mwtab_files/ST000122_AN000204_duplicate_subsections.txt", "r", encoding="utf-8") as f:
        mwtabfile.read(f)
    
    assert mwtabfile['METABOLOMICS WORKBENCH']['STUDY_ID'] == 'ST000123'
    assert mwtabfile['STUDY']['NUM_GROUPS'] == 'NA NA'
    


def test_write():
    mwtabfile = mwtab.mwtab.MWTabFile("tests/example_data/other_mwtab_files/ST000122_AN000204_duplicate_keys.txt")
    with open("tests/example_data/other_mwtab_files/ST000122_AN000204_duplicate_keys.txt", "r", encoding="utf-8") as f:
        mwtabfile.read(f)
    
    if not os.path.exists(os.path.dirname("tests/example_data/tmp/tmp.json")):
        dirname = os.path.dirname("tests/example_data/tmp/tmp.json")
        if dirname:
            os.makedirs(dirname)
    with open("tests/example_data/tmp/tmp.json", "w", encoding="utf-8") as f:
        mwtabfile.write(f, file_format="json")
    
    mwtabfile2 = mwtab.mwtab.MWTabFile("tests/example_data/tmp/tmp.json")
    with open("tests/example_data/tmp/tmp.json", "r", encoding="utf-8") as f:
        mwtabfile2.read(f)
    assert mwtabfile['MS_METABOLITE_DATA']['Data'][0]['Metabolite'] == '17-hydroxypregnenolone'
    
    
    mwtabfile = mwtab.mwtab.MWTabFile("tests/example_data/other_mwtab_files/ST000022_AN000041.txt")
    with open("tests/example_data/other_mwtab_files/ST000022_AN000041.txt", "r", encoding="utf-8") as f:
        mwtabfile.read(f)
    with open("tests/example_data/tmp/tmp.json", "w", encoding="utf-8") as f:
        mwtabfile.write(f, file_format="json")
    
    mwtabfile2 = mwtab.mwtab.MWTabFile("tests/example_data/tmp/tmp.json")
    with open("tests/example_data/tmp/tmp.json", "r", encoding="utf-8") as f:
        mwtabfile2.read(f)
    assert mwtabfile['NMR_BINNED_DATA']['Data'][0]['Metabolite'] == '0.4...0.46'
    
    
    with pytest.raises(IOError, match = r'^"filehandle" parameter must be writable'):
        with open("tests/example_data/tmp/tmp.json", "r", encoding="utf-8") as f:
            mwtabfile.write(f, file_format="json")
    
    with pytest.raises(TypeError, match = r'^Unknown file format'):
        with open("tests/example_data/tmp/tmp.json", "w", encoding="utf-8") as f:
            mwtabfile.write(f, file_format="asdf")
    
    with pytest.raises(TypeError, match = r'^Unknown file format'):
        mwtabfile.writestr(file_format="asdf")

    
def test_writing_results_file():
    results_file_dict = {'filename': 'ST000071_AN000111_Results.txt',
                          'UNITS': 'Peak area',
                          'Has m/z': 'Yes',
                          'Has RT': 'Yes',
                          'RT units': 'Minutes'}
    
    if not os.path.exists(os.path.dirname("tests/example_data/tmp/tmp.json")):
        dirname = os.path.dirname("tests/example_data/tmp/tmp.json")
        if dirname:
            os.makedirs(dirname)
    mwtabfile = mwtab.mwtab.MWTabFile("tests/example_data/other_mwtab_files/ST000022_AN000041_results_file.txt")
    with open("tests/example_data/other_mwtab_files/ST000022_AN000041_results_file.txt", "r", encoding="utf-8") as f:
        mwtabfile.read(f)
    
    # write json
    with open("tests/example_data/tmp/tmp.json", "w", encoding="utf-8") as f:
        mwtabfile.write(f, file_format="json")
    
    mwtabfile2 = mwtab.mwtab.MWTabFile("tests/example_data/tmp/tmp.json")
    with open("tests/example_data/tmp/tmp.json", "r", encoding="utf-8") as f:
        mwtabfile2.read(f)
    assert mwtabfile2['NM']['NMR_RESULTS_FILE'] == results_file_dict
    
    # write mwtab
    with open("tests/example_data/tmp/tmp.txt", "w", encoding="utf-8") as f:
        mwtabfile.write(f, file_format="mwtab")
    
    mwtabfile3 = mwtab.mwtab.MWTabFile("tests/example_data/tmp/tmp.txt")
    with open("tests/example_data/tmp/tmp.txt", "r", encoding="utf-8") as f:
        mwtabfile3.read(f)
    assert mwtabfile3['NM']['NMR_RESULTS_FILE'] == results_file_dict
    
    # write results file line missing filename
    del mwtabfile['NM']['NMR_RESULTS_FILE']['filename']
    del results_file_dict['filename']
    with open("tests/example_data/tmp/tmp.txt", "w", encoding="utf-8") as f:
        mwtabfile.write(f, file_format="mwtab")
    
    mwtabfile4 = mwtab.mwtab.MWTabFile("tests/example_data/tmp/tmp.txt")
    with open("tests/example_data/tmp/tmp.txt", "r", encoding="utf-8") as f:
        mwtabfile4.read(f)
    assert mwtabfile4['NM']['NMR_RESULTS_FILE'] == results_file_dict


def test_write_extended_section():
    mwtabfile = mwtab.mwtab.MWTabFile("tests/example_data/other_mwtab_files/ST000022_AN000041.txt")
    with open("tests/example_data/other_mwtab_files/ST000022_AN000041.txt", "r", encoding="utf-8") as f:
        mwtabfile.read(f)
    
    if not os.path.exists(os.path.dirname("tests/example_data/tmp/tmp.json")):
        dirname = os.path.dirname("tests/example_data/tmp/tmp.json")
        if dirname:
            os.makedirs(dirname)
    
    mwtabfile['NMR_BINNED_DATA']['Extended'] = []
    df = pandas.DataFrame([['some name', 'some value'], ['some name 2', 'some value 2']], columns = ['Metabolite', 'column1'])
    mwtabfile.set_extended_from_pandas(df, True)
    with open("tests/example_data/tmp/tmp.json", "w", encoding="utf-8") as f:
        mwtabfile.write(f, file_format="json")
    
    mwtabfile2 = mwtab.mwtab.MWTabFile("tests/example_data/tmp/tmp.json")
    with open("tests/example_data/tmp/tmp.json", "r", encoding="utf-8") as f:
        mwtabfile2.read(f)
    assert mwtabfile2['NMR_BINNED_DATA']['Extended'] == [{'Metabolite': 'some name', 'column1': 'some value'}, 
                                                          {'Metabolite': 'some name 2', 'column1': 'some value 2'}]
    
    
    with open("tests/example_data/tmp/tmp.txt", "w", encoding="utf-8") as f:
        mwtabfile.write(f, file_format="mwtab")
    
    mwtabfile3 = mwtab.mwtab.MWTabFile("tests/example_data/tmp/tmp.txt")
    with open("tests/example_data/tmp/tmp.txt", "r", encoding="utf-8") as f:
        mwtabfile3.read(f)
    assert mwtabfile3['NMR_BINNED_DATA']['Extended'] == [{'Metabolite': 'some name', 'column1': 'some value'}, 
                                                          {'Metabolite': 'some name 2', 'column1': 'some value 2'}]
    
    
    df = pandas.DataFrame([['some name', 'some value'], ['some name 3', 'some value 3']], columns = ['Metabolite', 'column2'])
    mwtabfile.set_extended_from_pandas(df)
    with open("tests/example_data/tmp/tmp.txt", "w", encoding="utf-8") as f:
        mwtabfile.write(f, file_format="mwtab")
    
    mwtabfile3 = mwtab.mwtab.MWTabFile("tests/example_data/tmp/tmp.txt")
    with open("tests/example_data/tmp/tmp.txt", "r", encoding="utf-8") as f:
        mwtabfile3.read(f)
    assert mwtabfile3['NMR_BINNED_DATA']['Extended'] == [{'Metabolite': 'some name', 'column2': 'some value'}, 
                                                          {'Metabolite': 'some name 3', 'column2': 'some value 3'}]
    
    
    mwtabfile['NMR_BINNED_DATA']['Extended'] = []
    mwtabfile._extended_metabolite_header = None
    with open("tests/example_data/tmp/tmp.txt", "w", encoding="utf-8") as f:
        mwtabfile.write(f, file_format="mwtab")
    
    mwtabfile3 = mwtab.mwtab.MWTabFile("tests/example_data/tmp/tmp.txt")
    with open("tests/example_data/tmp/tmp.txt", "r", encoding="utf-8") as f:
        mwtabfile3.read(f)
    assert mwtabfile3['NMR_BINNED_DATA']['Extended'] == []


def test_write_edge_cases(capsys):
    """Just hitting a few lines that are somewhat edge cases."""
    mwtabfile = mwtab.mwtab.MWTabFile("tests/example_data/other_mwtab_files/ST000022_AN000041.txt")
    with open("tests/example_data/other_mwtab_files/ST000022_AN000041.txt", "r", encoding="utf-8") as f:
        mwtabfile.read(f)
    
    if not os.path.exists(os.path.dirname("tests/example_data/tmp/tmp.json")):
        dirname = os.path.dirname("tests/example_data/tmp/tmp.json")
        if dirname:
            os.makedirs(dirname)
    
    save_binned_header = copy.deepcopy(mwtabfile._binned_header)
    mwtabfile._binned_header = None
    with open("tests/example_data/tmp/tmp.json", "w", encoding="utf-8") as f:
        mwtabfile.write(f, file_format="json")
    
    mwtabfile2 = mwtab.mwtab.MWTabFile("tests/example_data/tmp/tmp.json")
    with open("tests/example_data/tmp/tmp.json", "r", encoding="utf-8") as f:
        mwtabfile2.read(f)
    assert mwtabfile2._binned_header == save_binned_header
    
    
    with open("tests/example_data/tmp/tmp.txt", "w", encoding="utf-8") as f:
        mwtabfile.write(f, file_format="mwtab")
    
    mwtabfile3 = mwtab.mwtab.MWTabFile("tests/example_data/tmp/tmp.txt")
    with open("tests/example_data/tmp/tmp.txt", "r", encoding="utf-8") as f:
        mwtabfile3.read(f)
    assert mwtabfile3._binned_header == save_binned_header
    
    
    mwtabfile['NMR_BINNED_DATA']['Data'][0]['Bin range(ppm)'] = 'asdf'
    with open("tests/example_data/tmp/tmp.json", "w", encoding="utf-8") as f:
        mwtabfile.write(f, file_format="json")
    captured = capsys.readouterr()
    assert captured.out[:-1] == ("Warning: The \"Metabolite\" key and \"Bin range(ppm)\" "
                                  "key in ['NMR_BINNED_DATA']['Data'][0] are different "
                                  "values. Only the value in \"Bin range(ppm)\" will be written out.")
    
    
    mwtabfile4 = mwtab.mwtab.MWTabFile("tests/example_data/other_mwtab_files/ST000122_AN000204_extra_sample.txt")
    with open("tests/example_data/other_mwtab_files/ST000122_AN000204_extra_sample.txt", "r", encoding="utf-8") as f:
        mwtabfile4.read(f)
    assert mwtabfile4._factors is not None
        
    save_samples = copy.deepcopy(mwtabfile4._samples)
    mwtabfile4._samples = None
    mwtabfile4._factors = None
    with open("tests/example_data/tmp/tmp.json", "w", encoding="utf-8") as f:
        mwtabfile4.write(f, file_format="json")
    
    mwtabfile5 = mwtab.mwtab.MWTabFile("tests/example_data/tmp/tmp.json")
    with open("tests/example_data/tmp/tmp.json", "r", encoding="utf-8") as f:
        mwtabfile5.read(f)
    assert mwtabfile5._samples == save_samples
    assert mwtabfile5._factors is None
    
    
    with open("tests/example_data/tmp/tmp.txt", "w", encoding="utf-8") as f:
        mwtabfile4.write(f, file_format="mwtab")
    
    mwtabfile6 = mwtab.mwtab.MWTabFile("tests/example_data/tmp/tmp.txt")
    with open("tests/example_data/tmp/tmp.txt", "r", encoding="utf-8") as f:
        mwtabfile6.read(f)
    assert mwtabfile6._samples == save_samples
    assert mwtabfile6._factors is None

def test_write_error(init_tmp_dir):
    mwtabfile = mwtab.mwtab.MWTabFile("tests/example_data/other_mwtab_files/ST000122_AN000204_bad_table_type.json", force=True)
    with open("tests/example_data/other_mwtab_files/ST000122_AN000204_bad_table_type.json", "r", encoding="utf-8") as f:
        mwtabfile.read(f)
    mwtabfile['PROJECT'] = []
    with pytest.raises(TypeError, match = r'^Key/section "PROJECT" is not a dictionary\. It cannot be translated to the mwTab format\.'):
        with open("tests/example_data/tmp/tmp.txt", "w", encoding="utf-8") as f:
            mwtabfile.write(f, file_format="mwtab")
    

def test_keys_reorder():
    """Test that the keys are reordered to match what the Workbench wants."""
    test_dict = {"UNKNOWN_KEY": {'asdf': 'qwer'},
                  "MS_METABOLITE_DATA": {"Data":[{"attribute1":"qwer", "attribute2":"zxcv", "Metabolite":"asdf"}], "Units":"some unit"}, 
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
                                                    'attribute2': 'zxcv'}]},
                      "UNKNOWN_KEY": {'asdf': 'qwer'}}

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


@pytest.mark.parametrize("file_source", [
    "tests/example_data/other_mwtab_files/ST000122_AN000204_duplicate_keys.txt",
    "tests/example_data/other_mwtab_files/ST000122_AN000204_duplicate_keys.json"
])
def test_read_in_duplicate_keys(file_source):
    """Test that a file with duplicate keys is handled correctly."""
    
    if not os.path.exists("tests/example_data/tmp/"):
        os.makedirs("tests/example_data/tmp/")
    
    mwtabfile = mwtab.mwtab.MWTabFile(file_source, duplicate_keys=True)
    with open(file_source, "r", encoding="utf-8") as f:
        mwtabfile.read(f)
    
    assert 'key_1{{{_1_}}}' in mwtabfile["SUBJECT_SAMPLE_FACTORS"][0]["Additional sample data"]
    
    if file_source.endswith('.txt'):
        outpath = "tests/example_data/tmp/tmp.txt"
        file_format = 'mwtab'
    else:
        outpath = "tests/example_data/tmp/tmp.json"
        file_format = 'json'
    with open(outpath, "w", encoding="utf-8") as f:
        mwtabfile.write(f, file_format=file_format)
    
    new_mwtabfile = mwtab.mwtab.MWTabFile(outpath, duplicate_keys=True)
    with open(outpath, "r", encoding="utf-8") as f:
        new_mwtabfile.read(f)
    
    assert 'key_1{{{_1_}}}' in new_mwtabfile["SUBJECT_SAMPLE_FACTORS"][0]["Additional sample data"]



def test_validate():
    """Test that the validate method validates the object."""
    
    mwtabfile = mwtab.mwtab.MWTabFile("tests/example_data/other_mwtab_files/ST000122_AN000204_duplicate_keys.txt", duplicate_keys=True)
    with open("tests/example_data/other_mwtab_files/ST000122_AN000204_duplicate_keys.txt", "r", encoding="utf-8") as f:
        mwtabfile.read(f)
    
    error_log, _ = mwtabfile.validate(verbose=False)
    
    assert "duplicate keys" in error_log
    
    
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
    
    
    assert mwtabfile.study_id == "asdf"
    assert mwtabfile.analysis_id == "qwer"
    
    with pytest.raises(ValueError, match=r'^Header cannot be set because it is not of the form'):
        mwtabfile.header = "zxcv"
        
    mwtabfile.header = "#METABOLOMICS WORKBENCH STUDY_ID:ST000234 ANALYSIS_ID:AN000234 PROJECT_ID:PR000234"
    
    assert mwtabfile.header == "#METABOLOMICS WORKBENCH STUDY_ID:ST000234 ANALYSIS_ID:AN000234 PROJECT_ID:PR000234"
    
    with pytest.raises(TypeError, match=r'^The value for study_id must be a string.'):
        mwtabfile.study_id = 1
    
    mwtabfile['METABOLOMICS WORKBENCH'] = 'asdf'
    with pytest.raises(TypeError, match=r'^The "METABOLOMICS WORKBENCH" key is not a dictionary, so study_id cannot be set.'):
        mwtabfile.study_id = 'asdf'
    
    del mwtabfile['METABOLOMICS WORKBENCH']
    mwtabfile.study_id = 'opui'
    assert mwtabfile.study_id == 'opui'
    assert mwtabfile['METABOLOMICS WORKBENCH'] == {'STUDY_ID': 'opui'}
    
    del mwtabfile['METABOLOMICS WORKBENCH']
    mwtabfile.header = "#METABOLOMICS WORKBENCH STUDY_ID:ST000234 ANALYSIS_ID:AN000234 PROJECT_ID:PR000234"
    assert mwtabfile.header == "#METABOLOMICS WORKBENCH STUDY_ID:ST000234 ANALYSIS_ID:AN000234 PROJECT_ID:PR000234"
    assert mwtabfile['METABOLOMICS WORKBENCH'] == {'STUDY_ID': 'ST000234', 'ANALYSIS_ID': 'AN000234', 'PROJECT_ID': 'PR000234'}
    
    with pytest.raises(AttributeError):
        del mwtabfile.study_id


def test_get_and_set_table_from_pandas():
    """Test that METABOLITES DATA, METABOLITES, and EXTENDED can be set from a pandas dataframe and gotten as one."""
    
    mwtabfile = mwtab.mwtab.MWTabFile("tests/example_data/other_mwtab_files/ST000122_AN000204_duplicate_keys.txt")
        
    with open("tests/example_data/other_mwtab_files/ST000122_AN000204_duplicate_keys.txt", "r", encoding="utf-8") as f:
        mwtabfile.read(f)
    
    df = pandas.DataFrame([['some name', 'some value'], ['some name 2', 'some value 2']], columns = ['Metabolite', 'column1'])
    mwtabfile.set_table_from_pandas(df, 'Metabolites')
    assert mwtabfile['MS_METABOLITE_DATA']['Metabolites'] == [{'Metabolite': 'some name', 'column1': 'some value'}, 
                                                              {'Metabolite': 'some name 2', 'column1': 'some value 2'}]
    assert mwtabfile.get_table_as_pandas('Metabolites').equals(df)
    assert mwtabfile.get_metabolites_as_pandas().equals(df)
    assert mwtabfile._metabolite_header == ['column1']
    
    mwtabfile.set_table_from_pandas(df, 'Extended')
    assert mwtabfile['MS_METABOLITE_DATA']['Extended'] == [{'Metabolite': 'some name', 'column1': 'some value'}, 
                                                            {'Metabolite': 'some name 2', 'column1': 'some value 2'}]
    assert mwtabfile.get_table_as_pandas('Extended').equals(df)
    assert mwtabfile.get_extended_as_pandas().equals(df)
    assert mwtabfile._extended_metabolite_header == ['column1']
    
    mwtabfile.set_table_from_pandas(df, 'Data')
    assert mwtabfile['MS_METABOLITE_DATA']['Data'] == [{'Metabolite': 'some name', 'column1': 'some value'}, 
                                                        {'Metabolite': 'some name 2', 'column1': 'some value 2'}]
    assert mwtabfile.get_table_as_pandas('Data').equals(df)
    assert mwtabfile.get_metabolites_data_as_pandas().equals(df)
    assert mwtabfile._samples == ['column1']
    
    df2 = pandas.DataFrame([['some name3', 'some value3'], ['some name 4', 'some value 4']], columns = ['Metabolite', 'column2'])
    mwtabfile.set_metabolites_from_pandas(df2, True)
    assert mwtabfile['MS_METABOLITE_DATA']['Metabolites'] == [{'Metabolite': 'some name3', 'column2': 'some value3'}, 
                                                              {'Metabolite': 'some name 4', 'column2': 'some value 4'}]
    assert mwtabfile._metabolite_header == None
    
    mwtabfile.set_extended_from_pandas(df2, True)
    assert mwtabfile['MS_METABOLITE_DATA']['Extended'] == [{'Metabolite': 'some name3', 'column2': 'some value3'}, 
                                                            {'Metabolite': 'some name 4', 'column2': 'some value 4'}]
    assert mwtabfile._extended_metabolite_header == None
    
    mwtabfile.set_metabolites_data_from_pandas(df2, True)
    assert mwtabfile['MS_METABOLITE_DATA']['Data'] == [{'Metabolite': 'some name3', 'column2': 'some value3'}, 
                                                        {'Metabolite': 'some name 4', 'column2': 'some value 4'}]
    assert mwtabfile._samples == None
    
    del mwtabfile['MS_METABOLITE_DATA']
    mwtabfile['NMR_BINNED_DATA'] = {'Data':[{'Metabolite': 'asdf', 'sample1': '1234.45'}]}
    mwtabfile.set_metabolites_data_from_pandas(df2)
    assert mwtabfile['NMR_BINNED_DATA']['Data'] == [{'Metabolite': 'some name3', 'column2': 'some value3'}, 
                                                    {'Metabolite': 'some name 4', 'column2': 'some value 4'}]
    assert mwtabfile._samples == ['column2']
    assert mwtabfile._binned_header == ['column2']
    
    mwtabfile.set_metabolites_data_from_pandas(df, True)
    assert mwtabfile['NMR_BINNED_DATA']['Data'] == [{'Metabolite': 'some name', 'column1': 'some value'}, 
                                                    {'Metabolite': 'some name 2', 'column1': 'some value 2'}]
    assert mwtabfile._samples == None
    assert mwtabfile._binned_header == None
    

def test_print_file():
    mwtabfile = mwtab.mwtab.MWTabFile("tests/example_data/other_mwtab_files/ST000022_AN000041.txt")
    with open("tests/example_data/other_mwtab_files/ST000022_AN000041.txt", "r", encoding="utf-8") as f:
        mwtabfile.read(f)
    
    mwtab_io = io.StringIO()
    mwtabfile.print_file(mwtab_io, 'json')
    mwtab_str = mwtab_io.getvalue()
    assert mwtabfile.writestr('json') == mwtab_str[:-1]


def test_print_block():
    mwtabfile = mwtab.mwtab.MWTabFile("tests/example_data/other_mwtab_files/ST000022_AN000041.txt")
    with open("tests/example_data/other_mwtab_files/ST000022_AN000041.txt", "r", encoding="utf-8") as f:
        mwtabfile.read(f)
    
    mwtab_io = io.StringIO()
    mwtabfile.print_block('PROJECT', mwtab_io, 'json')
    mwtab_str = mwtab_io.getvalue()
    check_str = json.dumps(mwtabfile['PROJECT'], sort_keys=mwtab.mwtab.SORT_KEYS, indent=mwtab.mwtab.INDENT)
    assert check_str == mwtab_str[:-1]


def test_copy():
    """Copy and deepcopy dunders were overweritten, so check them."""
    mwtabfile = mwtab.mwtab.MWTabFile("tests/example_data/other_mwtab_files/ST000022_AN000041.txt")
    with open("tests/example_data/other_mwtab_files/ST000022_AN000041.txt", "r", encoding="utf-8") as f:
        mwtabfile.read(f)
    
    mwtabfile2 = copy.copy(mwtabfile)
    assert mwtabfile == mwtabfile2
    assert mwtabfile._samples == mwtabfile2._samples
    assert mwtabfile._duplicate_keys == mwtabfile2._duplicate_keys
    
    mwtabfile3 = copy.deepcopy(mwtabfile)
    assert mwtabfile == mwtabfile3
    assert mwtabfile._samples == mwtabfile3._samples
    assert mwtabfile._duplicate_keys == mwtabfile3._duplicate_keys


def test_misc_coverage():
    """Testing some lines that are basic error checking that don't really belong in another test."""
    mwtabfile = mwtab.mwtab.MWTabFile("tests/example_data/other_mwtab_files/ST000022_AN000041.txt")
    with open("tests/example_data/other_mwtab_files/ST000022_AN000041.txt", "r", encoding="utf-8") as f:
        mwtabfile.read(f)
    
    with pytest.raises(TypeError, match = r"^Expecting <class 'str'> or <class 'bytes'>, but <class 'dict'> was passed"):
        mwtabfile._is_mwtab({})
    
    with pytest.raises(TypeError, match = r"^Expecting <class 'str'> or <class 'bytes'>, but <class 'dict'> was passed"):
        mwtabfile._is_json({})



