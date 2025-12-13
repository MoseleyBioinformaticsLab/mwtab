import pytest
import mwtab
import copy


@pytest.mark.parametrize("file_source", [
    "tests/example_data/validation_files/ST000122_AN000204_validate_passing.txt",
    "tests/example_data/validation_files/ST000122_AN000204_validate_passing.json"
])
def test_validate_passing(file_source):
    mwfile = next(mwtab.read_files(file_source))
    validation_log, _ = mwtab.validate_file(mwfile)
    assert len(validation_log.split('\n')) == 9
    assert 'Status: Passing' in validation_log


@pytest.mark.parametrize("file_source", [
    "tests/example_data/validation_files/ST000122_AN000204_validate_polarity.txt",
    "tests/example_data/validation_files/ST000122_AN000204_validate_polarity.json"
])
def test_validate_polarity(file_source):
    mwfile = next(mwtab.read_files(file_source))
    validation_log, _ = mwtab.validate_file(mwfile)
    assert 'Error: The "polarity" column in the ' in validation_log
    assert ('indicates multiple polarities in a single analysis, and '
            'this should not be. A single mwTab file is supposed to be '
            'restricted to a single analysis. This means multiple MS '
            'runs under different settings should each be in their own file.') in validation_log
    

@pytest.mark.parametrize("file_source", [
    "tests/example_data/validation_files/ST000122_AN000204_validate_table_values.txt",
    "tests/example_data/validation_files/ST000122_AN000204_validate_table_values.json"
])
def test_validate_table_values(file_source):
    mwfile = next(mwtab.read_files(file_source))
    validation_log, _ = mwtab.validate_file(mwfile)
    assert 'Error: Column(s) with no name were found in the' in validation_log
    
    assert 'Warning: The "null_column" column at position' in validation_log
    assert 'table has all null values.' in validation_log
    
    assert 'Warning: The "90_10" column at position' in validation_log
    assert ('table may have incorrect values. 90% or more of the values are the '
            'same, but 10% or less are different.') in validation_log
    
    assert 'Warning: There are duplicate rows in the' in validation_log
    
    assert 'Warning: There are duplicate column names in the' in validation_log

def test_validate_table_values2():
    """Can only be tested for JSON files."""
    mwfile = next(mwtab.read_files("tests/example_data/validation_files/ST000122_AN000204_validate_table_values2.json"))
    validation_log, _ = mwtab.validate_file(mwfile)
    assert ' table does not have the same columns for every row.' in validation_log



@pytest.mark.parametrize("file_source", [
    "tests/example_data/validation_files/ST000122_AN000204_validate_metabolite_names.txt",
    "tests/example_data/validation_files/ST000122_AN000204_validate_metabolite_names.json"
])
def test_validate_metabolite_names(file_source):
    mwfile = next(mwtab.read_files(file_source))
    validation_log, _ = mwtab.validate_file(mwfile)
    # This one should bein the Data table.
    assert 'Warning: There is a metabolite name, "samples"' in validation_log
    # This one should be in the Metabolites table.
    assert 'Warning: There is a metabolite name, "factors"' in validation_log
    # This one should be in the Extended table.
    assert 'Warning: There is a metabolite name, "metabolite_name"' in validation_log
    
    # All 3 will have this message, but we only need to confirm this part once.
    assert (' that is probably wrong. It is close to a header name and '
            'is likely due to a badly constructed Tab file.') in validation_log



@pytest.mark.parametrize("file_source", [
    "tests/example_data/validation_files/ST000122_AN000204_validate_extended.txt",
    "tests/example_data/validation_files/ST000122_AN000204_validate_extended.json"
])
def test_validate_extended(file_source):
    mwfile = next(mwtab.read_files(file_source))
    validation_log, _ = mwtab.validate_file(mwfile)
    assert ' table has Sample IDs that were not found' in validation_log
    assert 'Those IDs are:\n\t"16_A0_Lung_naive_0days_170427_UKy_GCH_rep1-lipid"' in validation_log


@pytest.mark.parametrize("file_source", [
    "tests/example_data/validation_files/ST000122_AN000204_validate_extended2.txt",
    "tests/example_data/validation_files/ST000122_AN000204_validate_extended2.json"
])
def test_validate_extended2(file_source):
    """Cannot be tested at the same time as the other test for extended."""
    mwfile = next(mwtab.read_files(file_source))
    validation_log, _ = mwtab.validate_file(mwfile)
    assert ' table does not have a column for "sample_id".' in validation_log




@pytest.mark.parametrize("file_source", [
    "tests/example_data/validation_files/ST000122_AN000204_validate_ssf.txt",
    "tests/example_data/validation_files/ST000122_AN000204_validate_ssf.json"
])
def test_validate_subject_samples_factors(file_source):
    mwfile = next(mwtab.read_files(file_source))
    validation_log, _ = mwtab.validate_file(mwfile)
    assert ' has a duplicate Sample ID.' in validation_log
    
    assert 'has the following duplicate keys in its Factors:\n\t"Tissue/Fluid"' in validation_log
    
    assert 'has the following duplicate keys in its Additional sample data:\n\t"key1"' in validation_log


def test_validate_factors():
    """This can only be detected for mwTab, not JSON."""
    mwfile = next(mwtab.read_files("tests/example_data/validation_files/ST000122_AN000204_validate_factors.txt"))
    validation_log, _ = mwtab.validate_file(mwfile)
    assert ("Error: The factors in the METABOLITE_DATA section "
            "and SUBJECT_SAMPLE_FACTORS section do not match.") in validation_log





def test_validate_header_lengths():
    """This can only be detected for mwTab, not JSON."""
    mwfile = next(mwtab.read_files("tests/example_data/validation_files/ST000122_AN000204_validate_header_lengths.txt"))
    validation_log, _ = mwtab.validate_file(mwfile)
    assert ('Error: The section, METABOLITES, has a mismatch between the '
            'number of headers and the number of elements in each line. Either '
            'a line(s) has more values than headers or there are too few headers.') in validation_log
    assert ('Error: The section, MS_METABOLITE_DATA, has a mismatch between '
            'the number of headers and the number of elements in each line. '
            'Either a line(s) has more values than headers or there are too few headers.') in validation_log


def test_validate_sub_section_uniqueness():
    """This can only be detected for mwTab, not JSON."""
    mwfile = next(mwtab.read_files("tests/example_data/validation_files/ST000122_AN000204_validate_subsection_uniqueness.txt"))
    validation_log, _ = mwtab.validate_file(mwfile)
    assert "Error: The section, STUDY, has a sub-section, FIRST_NAME, that is duplicated." in validation_log


@pytest.mark.parametrize("file_source", [
    "tests/example_data/validation_files/ST000122_AN000204_validate_metabolites.txt",
    "tests/example_data/validation_files/ST000122_AN000204_validate_metabolites.json"
])
def test_validate_metabolites(file_source):
    mwfile = next(mwtab.read_files(file_source))
    validation_log, _ = mwtab.validate_file(mwfile)
    assert 'Warning: The following metabolites in the' in validation_log
    assert 'were not found in the' in validation_log
    assert 'table:\n\t"Corticosterone, DOC"' in validation_log
    
    assert 'Error: A metabolite without a name was found in the' in validation_log
    
    assert 'Warning: The following metabolites in the' in validation_log
    assert 'table appear more than once in the table:\n\t"Testosterone"' in validation_log
    
    assert 'Warning: The "ri" column at position ' in validation_log
    assert ' table, matches a standard column name, "retention_index"' in validation_log
    
    assert 'Warning: The "pubchem_id" column at position ' in validation_log
    assert 'and some of the values in the column do not match the expected type or format for that column.' in validation_log
    
    assert 'Warning: The column "ri" was found in the' in validation_log
    assert 'but this column implies that another column, "retention_index_type",' in validation_log
    
    assert 'Warning: The column pair, "other_id" and "other_id_type",' in validation_log
    assert ('in the METABOLITES table should have data in the '
            'same rows, but at least one row has data in one '
            'column and nothing in the other.') in validation_log
    
    assert 'Warning: The standard column, "other_id", was' in validation_log
    assert ('If this column contains database IDs for standard databases such '
            'as KEGG, PubChem, HMDB, etc., it is recommended to make individual '
            'columns for these and not lump them together into a less descriptive '
            '"other_id" column.') in validation_log
    
    assert 'Warning: The column, "pubchem_id/hmdb_id", in the' in validation_log
    assert ('was matched to multiple standard names, [\'pubchem_id\', \'hmdb_id\']. This is a good indication '
            'that the values in that column should be split into the appropriate '
            'individual columns.') in validation_log
    
        

@pytest.mark.parametrize("file_source", [
    "tests/example_data/validation_files/ST000122_AN000204_validate_metabolite_data.txt",
    "tests/example_data/validation_files/ST000122_AN000204_validate_metabolite_data.json"
])
def test_validate_data(file_source):
    mwfile = next(mwtab.read_files(file_source))
    validation_log, _ = mwtab.validate_file(mwfile)
    assert 'Error: SUBJECT_SAMPLE_FACTORS section missing sample ID(s).' in validation_log
    assert 'The following IDs were found in the' in validation_log 
    assert 'section but not in the SUBJECT_SAMPLE_FACTORS:\n\t"CER040_242995_ML_02"' in validation_log
    
    assert 'Warning: There are duplicate samples in the' in validation_log
    
    assert 'Warning: The following metabolites in the' in validation_log
    assert 'were not found in the' in validation_log
    assert 'table:\n\t"Corticosterone_ DOC"' in validation_log
    
    assert 'Error: A metabolite without a name was found in the' in validation_log
    
    assert 'Warning: The following metabolites in the' in validation_log
    assert 'table appear more than once in the table:\n\t"Testosterone"' in validation_log
    


@pytest.mark.parametrize("file_source", [
    "mwtab",
    "json"
])
class TestMWTabSchemaErrors:
    validation_logs = {}
    mwfile = next(mwtab.read_files("tests/example_data/validation_files/ST000122_AN000204_schema_errors.txt"))
    validation_log, _ = mwtab.validate_file(mwfile)
    validation_logs['mwtab'] = validation_log
    mwfile = next(mwtab.read_files("tests/example_data/validation_files/ST000122_AN000204_schema_errors.json"))
    validation_log, _ = mwtab.validate_file(mwfile)
    validation_logs['json'] = validation_log
    
    
    def test_ionization(self, file_source):
        """IONIZATION should not have a value of positive or negative, there is a special check for this."""
        assert '"ION_MODE" is where that should be indicated.' in self.validation_logs[file_source]
    
    def test_custom_message(self, file_source):
        """Some schemas can have custom messages, such as COLLISON_GAS."""
        assert 'should be one of "Nitrogen" or "Argon".' in self.validation_logs[file_source]
    
    def test_message(self, file_source):
        """Some schemas can have the entire message inserted. This is different from custom message."""
        assert ('Error: There must be either a "MS_METABOLITE_DATA" '
                'section or a "MS_RESULTS_FILE" subsection in the '
                '"MS" section. Neither were found.') in self.validation_logs[file_source]
    
    def test_additionalProperties_errors(self, file_source):
        """There are a few different places an additionalProperties can occur."""
        assert 'Unknown or invalid section, "BADSECTION".' in self.validation_logs[file_source]
        assert 'Unknown or invalid subsections, "BADSUBSECTION", "BADSUBSECTION2", ' in self.validation_logs[file_source]
        assert 'Unknown or invalid subsection, "BADSUBSECTION", ' in self.validation_logs[file_source]
    
    def test_required_error(self, file_source):
        """Simply test that an error is printed wwhen COLLECTION is not present."""
        assert 'The required property, "COLLECTION", ' in self.validation_logs[file_source]
        assert ' is missing.' in self.validation_logs[file_source]
    
    def test_email_format_error(self, file_source):
        """Test that badly formatted emails are caught."""
        assert ' is not a valid email.' in self.validation_logs[file_source]
    
    def test_NA_value_error(self, file_source):
        """Test that messages are printed about NA values where they should not be."""
        assert 'An empty value or a null value was detected' in self.validation_logs[file_source]
        # required properties have a slightly different message than non-required ones.
        assert 'A legitimate value should be provided for this required' in self.validation_logs[file_source]
        assert 'Either a legitimate value should be provided for this' in self.validation_logs[file_source]



class TestMWTabSchemaErrorsRare:
    """None of what is tested here is expected to show up during normal operation. 
    There are more custom validations in create_better_error_messages than what are 
    acutally in the schema used to validate an mwTab file. This test is to cover 
    those lines.
    """
    ms_schema = copy.deepcopy(mwtab.mwschema.ms_required_schema)
    ms_schema['properties']['FOR_TEST'] = {'type':'object', 'properties':{}}
    ms_schema['properties']['FOR_TEST']['properties']['MIN_PROPERTIES'] = {'type':'object', 'minProperties':1}
    del ms_schema['properties']['CHROMATOGRAPHY']['properties']['COLUMN_NAME']['not']
    ms_schema['properties']['CHROMATOGRAPHY']['properties']['COLUMN_NAME']['minLength'] = 1
    ms_schema['properties']['CHROMATOGRAPHY']['properties']['FLOW_GRADIENT']['minLength'] = 1000
    ms_schema['properties']['CHROMATOGRAPHY']['properties']['FLOW_RATE']['maxLength'] = 3
    ms_schema['properties']['FOR_TEST']['properties']['MIN_ITEMS'] = {'type':'array', 'minItems':1}
    ms_schema['properties']['FOR_TEST']['properties']['MIN_ITEMS2'] = {'type':'array', 'minItems':2}
    ms_schema['properties']['CHROMATOGRAPHY']['properties']['CHROMATOGRAPHY_SUMMARY']['type'] = ['string', 'null']
    ms_schema['properties']['FOR_TEST']['properties']['ENUM'] = {'type':'string', 'enum':['asdf', 'qwer']}
    ms_schema['properties']['FOR_TEST']['properties']['PATTERN'] = {'type':'string', 'pattern':r'^asdf$'}
    ms_schema['properties']['FOR_TEST']['properties']['MIN'] = {'type':'integer', 'minimum':15}
    ms_schema['properties']['FOR_TEST']['properties']['MAX'] = {'type':'integer', 'maximum':5}
    ms_schema['properties']['FOR_TEST']['properties']['UNIQUE_ITEMS'] = {'type':'array', 'uniqueItems':True}
    ms_schema['properties']['FOR_TEST']['dependentRequired'] = {'MIN':['asdf']}
    ms_schema['properties']['FOR_TEST']['properties']['NOT_ONEOF'] = {'type':'string', 'not':{'oneOf':[{'enum':['qwer']}, {'pattern':r'asdf'}]}}
    
    mwfile = next(mwtab.read_files("tests/example_data/validation_files/ST000122_AN000204_schema_errors.txt"))
    mwfile['FOR_TEST'] = {}
    mwfile['FOR_TEST']['MIN_PROPERTIES'] = {}
    mwfile['CHROMATOGRAPHY']['CHROMATOGRAPHY_SUMMARY'] = 1
    mwfile['CHROMATOGRAPHY']['CHROMATOGRAPHY_TYPE'] = 1
    mwfile['CHROMATOGRAPHY']['COLUMN_NAME'] = ''
    mwfile['FOR_TEST']['MIN_ITEMS'] = []
    mwfile['FOR_TEST']['MIN_ITEMS2'] = ['a']
    mwfile['FOR_TEST']['ENUM'] = 'zxcv'
    mwfile['FOR_TEST']['PATTERN'] = 'zxcv'
    mwfile['FOR_TEST']['MIN'] = 10
    mwfile['FOR_TEST']['MAX'] = 10
    mwfile['FOR_TEST']['UNIQUE_ITEMS'] = ['a', 'a']
    mwfile['FOR_TEST']['NOT_ONEOF'] = 'asdf'
    validation_log, _ = mwtab.validate_file(mwfile, ms_schema=ms_schema)
    
    
    def test_minProperties_error(self):
        assert ' "MIN_PROPERTIES", in the "FOR_TEST" section cannot be empty' in self.validation_log
        
    def test_minLength_error(self):
        assert 'cannot be an empty string.' in self.validation_log
        assert 'is too short' in self.validation_log
    
    def test_maxLength_error(self):
        assert 'is too long' in self.validation_log
    
    def test_minItems_error(self):
        assert ' "MIN_ITEMS", in the "FOR_TEST" section cannot be empty' in self.validation_log
        assert 'must have at least 2 items' in self.validation_log
    
    def test_type_error(self):
        assert 'is not any of the allowed types:' in self.validation_log
        assert 'is not of type' in self.validation_log
    
    def test_enum_error(self):
        assert '"ENUM", in the "FOR_TEST" section is not one of [\'asdf\', \'qwer\'].' in self.validation_log
    
    def test_pattern_error(self):
        assert 'does not match the regular expression pattern' in self.validation_log
    
    def test_minimum_error(self):
        assert 'must be greater than or equal to' in self.validation_log
    
    def test_maximum_error(self):
        assert 'must be less than or equal to' in self.validation_log
    
    def test_uniqueItems_error(self):
        assert 'has non-unique elements.' in self.validation_log
    
    def test_not_oneOf_error(self):
        assert '"NOT_ONEOF"' in self.validation_log
        assert ("does not match the regular expression pattern "
                "{'oneOf': [{'enum': ['qwer']}, {'pattern': 'asdf'}]}") in self.validation_log



@pytest.mark.parametrize("file_source", [
    "tests/example_data/mwtab_files/ST000122_AN000204.json"
])
def test_validation_log_local(file_source):
    mwfile = next(mwtab.read_files(file_source))
    validation_log, _ = mwtab.validate_file(mwfile)
    # assert "mwtab version: {}".format(mwtab.__version__) in validation_log
    assert "Source:        {}".format(file_source) in validation_log
    assert "Study ID:      {}".format("ST000122") in validation_log
    assert "Analysis ID:   {}".format("AN000204") in validation_log
    assert "File format:   {}".format("json") in validation_log


@pytest.mark.parametrize("file_source", [
    "2"
])
def test_validation_log_web(file_source):
    mwfile = next(mwtab.read_files(file_source))
    validation_log, _ = mwtab.validate_file(mwfile)
    # assert "mwtab version: {}".format(mwtab.__version__) in validation_log
    assert "Source:        {}".format("https://www.metabolomicsworkbench.org/rest/study/analysis_id/AN000002/mwtab/txt")\
            in validation_log
    assert "Study ID:      {}".format("ST000002") in validation_log
    assert "Analysis ID:   {}".format("AN000002") in validation_log
    assert "File format:   {}".format("txt") in validation_log


def test_base_schema_error_extra_key():
    """Test that if extra keys are included there is an error."""
    mwfile = next(mwtab.read_files("tests/example_data/mwtab_files/ST000122_AN000204.json"))
    mwfile["NM"] = {}
    validation_log, _ = mwtab.validate_file(mwfile)
        
    assert 'Error: Unknown or invalid sections, "CHROMATOGRAPHY", "MS", "MS_METABOLITE_DATA".' in validation_log
    assert ('Error: There must be either a "NMR_METABOLITE_DATA" section, a '
            '"NMR_BINNED_DATA" section or a "NMR_RESULTS_FILE" subsection in the "NM" section. Neither were found.') in validation_log
    

@pytest.mark.parametrize("file_source", [
    "tests/example_data/validation_files/ST000122_AN000204_validate_file.txt",
    "tests/example_data/validation_files/ST000122_AN000204_validate_file.json"
])
def test_validate_file(file_source):
    mwfile = next(mwtab.read_files(file_source))
    validation_log, _ = mwtab.validate_file(mwfile)
    assert ('Error: No "MS" or "NM" section was found, '
            'so analysis type could not be determined. '
            'Mass spec will be assumed.') in validation_log
    
    if file_source.endswith('.txt'):
        assert 'Warning: Missing METABOLITES section.' in validation_log
    else:
        assert 'Warning: Missing ["MS_METABOLITE_DATA"]["Metabolites"] section.' in validation_log
    

def test_validation_complete_coverage():
    """This is just to hit some lines that aren't covered, but also aren't terribly important to test."""
    mwfile = next(mwtab.read_files("tests/example_data/validation_files/complete_coverage.json"))
    validation_log, _ = mwtab.validate_file(mwfile, verbose=True)
    assert validation_log is None

def test_validation_complete_coverage2():
    """This is just to hit some lines that aren't covered, but also aren't terribly important to test."""
    mwfile = next(mwtab.read_files("tests/example_data/validation_files/complete_coverage2.txt"))
    validation_log, _ = mwtab.validate_file(mwfile)

def test_validation_complete_coverage3():
    """This is just to hit some lines that aren't covered, but also aren't terribly important to test."""
    mwfile = next(mwtab.read_files("tests/example_data/validation_files/complete_coverage3.json"))
    validation_log, _ = mwtab.validate_file(mwfile)
