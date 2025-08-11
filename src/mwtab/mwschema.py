#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
mwtab.mwschema
~~~~~~~~~~~~~~

This module provides schema definitions for different sections of the
``mwTab`` Metabolomics Workbench format.
"""

import sys
import re
from copy import deepcopy
from collections import OrderedDict
from functools import partial, partialmethod
from collections.abc import Callable

from schema import Schema, Optional, Or, And, SchemaError
import jsonschema

from . import metadata_column_matching

class _duplicate_key_list(list):
    """Class identical to list that can be used for type checking. Used to handle dealing with parsing duplicate keys in JSON."""
    def __init__(self, *args, **kwargs):
        super(_duplicate_key_list, self).__init__(*args, **kwargs)



if sys.version_info.major == 2:
    str = unicode


def create_units_regex(units: list[str], can_be_range: bool = False) -> str:
    """Create a regular expression to match something like '5 V' or '5-6 V'.
    
    Regular expression will match something like '5 V'. The space is required, '5V' 
    will not match. If can_be_range is True, then the expression will also match 
    something like '5-7 V' (hyphen), '5−7 V' (em dash), and '5 to 7 V' as well. 
    The spacing must match.
    
    Args:
        units: The list of unit strings that could follow the number.
        can_be_range: If True, allow a pattern like '5-7 V' with a range.
    
    Returns:
        A regular expression to match a number and units based on the arguments.
    """
    if can_be_range:
        regex = '^(' + metadata_column_matching.NUM_RANGE + '|' + metadata_column_matching.NUMS + ')' + f' ({"|".join(units)})$'
    else:
        regex = '^' + metadata_column_matching.NUMS + f' ({"|".join(units)})$'
    return regex

def create_unit_error_message(validation_error: jsonschema.exceptions.ValidationError, can_be_range: bool = False, no_units: bool = False, units: list[str]|None = None) -> str:
    """Generate the error message for mwTab subsections that fail the unit regex.
    
    The idea for this function is that you include it in the jsonschema under a keyword that 
    is not reserved by jsonschema, and then your custom error handler will look for that 
    keyword and execute whatever function it finds there, passing the error object in as 
    the first parameter. You can combine this with the "partial" function from functools 
    to set the other parameters of this function as needed.
    
    Args:
        validation_error: The ValidationError created by jsonschema. Used to fill in some text in the message.
        can_be_range: If True, the regular expression used to validate could have matched a number range, so the message is modified to note that.
        no_units: If True, the regular expression did not require units to be present, so the message is modified to note that.
        units: If the regular expression required units, pass them in with this parameter so the message will indicate the allowed units.
    
    Returns:
        A completed string error message.
    """
    section = validation_error.relative_path[-2]
    subsection = validation_error.relative_path[-1]
    
    if can_be_range:
        range_string = ' or range (ex. 5-6) '
    else:
        range_string = ''
    
    if no_units:
        unit_string = '.'
    else:
        unit_string = f'followed by a space with a unit from the following list: {units}.'
    
    message = (f'The value, "{validation_error.instance}", for the subsection, "{subsection}", in the "{section}" section '
               f'should be a number{range_string}{unit_string}')
    return message

def _create_unit_regex_and_message_func(units: list[str], can_be_range: bool = False) -> str:
    """Simple wrapper for DRY purposes.
    
    Creating the regular expression and validation error message in 1 function mixes too many 
    concerns into 1 function, so they are split into 2 and this function serves as a convenience 
    to pass them into a jsonschema easily. To this end the return is in dictionary form with the 
    intention for it to be unpacked into the jsonschema. For example, 
    {**_create_unit_regex_and_message_func(['V'], True)}
    
    Args:
        units: A list of strings used to create the regex and error message.
        can_be_range: If true, the regex will match a number range and the message will be slightly different.
    
    Returns:
        A dicitonary {'pattern': regex, 'message_func': message_function}.
    """
    regex = create_units_regex(units, can_be_range)
    message_function = partial(create_unit_error_message, can_be_range = can_be_range, no_units = False, units = units)
    return {'pattern': regex, 'message_func': message_function}








def create_prefixed_unit_list(unit: str, prefixes: list[str] = ['m', 'k', 'n', 'u', 'c']) -> list[str]:
    """Adds given prefixes to the given unit and returns a list of the concatentaions.
    
    Args:
        unit: The unit to add prefixes to.
        prefixes: The prefixes to concat to unit.
    
    Returns:
        A list of units prefixed with prefixes and the unit unprefixed.
    """
    concatenations = [prefix + unit for prefix in prefixes]
    return concatenations + [unit]
    

class UnitChecker(object):
    """A Validatable created to use with the Schema package.
    
    Validates whether a data element is a number followed by a space, then the given units or not.
    Recommended to be used like And(str, UnitChecker('Hz')).
    
    Parameters:
        units: The list of unit strings that could follow the number.
        can_be_range: If True, allow a pattern like '5-7 V' with a range.
        error_message: A custom error message to use instead of the default one. 
          This string will have the format method called on it with locals(), so variables such as 
          'self.units' and 'data' can be used in the message.
    
    Attributes:
        units (list[str]): The current list of unit strings that must follow the number.
        can_be_range (bool): Whether the checker allows a range (ex. '5-7 V') or not.
        error_message (str|None): The current custom error message to use instead of the default one.
    """
    def __init__(self, units: list[str], can_be_range: bool = False, 
                 error_message: str|None = None, error_message_builder: Callable|None = None, **kwargs):
        self.units = units
        self.can_be_range = can_be_range
        self.error_message = error_message
        self.error_message_builder = error_message_builder
        self.kwargs = kwargs

    def validate(self, data: str) -> bool:
        """Test if data is a number followed by a unit.
        
        data must be something like '5 V' to return True. The space is required, '5V' 
        will raise a SchemaError. If the can_be_range attribute is True, then data 
        can look like '5-7 V' as well.
        
        Args:
            data: the value to test, should be a string.
        
        Returns:
            True if the data is a unit, otherwise raises a SchemaError.
        """
        if isinstance(data, str):
            if self.can_be_range:
                regex = '(' + metadata_column_matching.NUM_RANGE + '|' + metadata_column_matching.NUMS + ')' + f' ({"|".join(self.units)})'
            else:
                regex = metadata_column_matching.NUMS + f' ({"|".join(self.units)})'
            if re.fullmatch(regex, data.strip()):
                return True
        if self.error_message_builder is not None:
            error_message = self.error_message_builder(self, data)
        elif self.error_message:
            error_message = self.error_message
        else:
            error_message = f'"{data}" should be a number followed by a space with a unit from the following list: {self.units}.'
        raise SchemaError(error_message.format(**locals()))

class NumChecker(object):
    """A Validatable created to use with the Schema package.
    
    Validates whether a data element is a number.
    Recommended to be used like And(str, NumChecker()).
    
    Parameters:
        error_message: A custom error message to use instead of the default one. 
          This string will have the format method called on it with locals(), so variables such as 
          'data' can be used in the message.
    
    Attributes:
        error_message (str|None): The current custom error message to use instead of the default one.
    """
    def __init__(self, error_message: str|None = None):
        self.error_message = error_message

    def validate(self, data: str) -> bool:
        """Test if data is a number.
        
        data must be something like '5' to return True.
        
        Args:
            data: the value to test, should be a string.
        
        Returns:
            True if the data is a number, otherwise raises a SchemaError.
        """
        if isinstance(data, str):
            if re.fullmatch(r'\d+', data.strip()):
                return True
        if self.error_message:
            error_message = self.error_message
        else:
            error_message = f'"{data}" should be a number.'
        raise SchemaError(error_message.format(**locals()))

NA_VALUES = ['', '-', '−', '--', '---', 
             'NA', 'na', 'n.a.', 'N.A.', 'n/a', 'N/A', '#N/A', 'NaN', 'nan', 
             'null', 'Null', 'NULL', 'none', 'None',
             'unspecified', 'Unspecified']
class SetChecker(object):
    f"""A Validatable created to use with the Schema package.
    
    Validates whether a data element is from a set of values.
    Recommended to be used like And(str, SetChecker()).
    
    Parameters:
        in_values: The set of values that the data element must be in to be valid.
        na_values: The set of NA values to look for and disallow. If set to None, 
          Defaults to a set of values pulled from Metabolomics Workbench datasets. 
          To ignore NA values, set this to an empty list.
          The set of default NA values are: 
              {NA_VALUES}
        error_message: A custom error message to use instead of the default one. 
          This string will have the format method called on it with locals(), so variables such as 
          'data' can be used in the message.
    
    Attributes:
        in_values (list[str]): The current set of values that a data element must be in to be valid.
        na_values (list[str]): The current set of NA values to disallow.
        error_message (str|None): The current custom error message to use instead of the default one.
    """
    def __init__(self, in_values: list[str], na_values: list[str]|None = None, 
                 error_message: str|None = None, error_message_builder: Callable|None = None, **kwargs):
        self.in_values = in_values
        self.na_values = NA_VALUES if na_values is None else na_values
        self.error_message = error_message
        self.error_message_builder = error_message_builder
        self.kwargs = kwargs

    def validate(self, data: str) -> bool:
        """Test if data is in in_values.
        
        Args:
            data: the value to test, should be a string.
        
        Returns:
            True if the data is in in_values, otherwise raises a SchemaError.
        """
        if isinstance(data, str):
            stripped_data = data.strip()
            if stripped_data in self.in_values:
                return True
        # if self.error_message:
        #     error_message = self.error_message
        # else:
        #     error_message = f'"{data}" should be one of the following values: {self.in_values}.'
        if self.error_message_builder is not None:
            error_message = self.error_message_builder(self, data)
        else:
            error_message = self._build_error_message(data)
        raise SchemaError(error_message.format(**locals()))
    
    def _build_error_message(self, data):
        """
        """
        if self.error_message:
            error_message = self.error_message
        else:
            error_message = f'"{data}" should be one of the following values: {self.in_values}.'
        return error_message

# def is_unit(value: str, unit: str):
#     """Test if value is a number followed by a unit.
    
#     value must be something like '5 V' to return True. The space is required, '5V' 
#     will return False.
    
#     Args:
#         value: the value to test, should be a string.
#         unit: the required unit string that should follow the number.
    
#     Returns:
#         True if the value is a unit, False otherwise.
#     """
#     if isinstance(value, str):
#         if re.match(metadata_column_matching.NUMS + ' ' + unit, value.strip()):
#             return True
#     return False

def create_error_message(section: str, subsection: str, can_be_range: bool = False, no_units: bool = False) -> str:
    """
    """
    if can_be_range:
        range_string = ' or range (ex. 5-6)'
    else:
        range_string = ''
    
    if no_units:
        unit_string = '.'
    else:
        unit_string = ' followed by a space with a unit from the following list: {self.units}.'
    
    message = ('The value, "{data}", for the subsection, ' + f'"{subsection}", in the "{section}" section '
               f'should be a number{range_string}{unit_string}')
    return message

def create_error_message2(self, data) -> str:
    """
    """
    can_be_range = self.can_be_range
    no_units = self.kwargs['no_units'] if 'no_units' in self.kwargs else False
    section = self.kwargs['section']
    subsection = self.kwargs['subsection']
    
    if can_be_range:
        range_string = ' or range (ex. 5-6) '
    else:
        range_string = ''
    
    if no_units:
        unit_string = '.'
    else:
        unit_string = 'followed by a space with a unit from the following list: {self.units}.'
    
    message = (f'The value, "{data}", for the subsection, "{subsection}", in the "{section}" section '
               f'should be a number{range_string}{unit_string}')
    return message

def partialclass(cls, *args, **kwargs):

    class NewCls(cls):
        __init__ = partialmethod(cls.__init__, *args, **kwargs)
        
        # def _build_error_message(self, data, **kwargs):
        #     return create_error_message2(self, data)
            
            
        #     can_be_range = self.can_be_range
        #     no_units = kwargs['no_units'] if 'no_units' in kwargs else False
        #     section = kwargs['section']
        #     subsection = kwargs['subsection']
            
        #     if can_be_range:
        #         range_string = ' or range (ex. 5 - 6) '
        #     else:
        #         range_string = ''
            
        #     if no_units:
        #         unit_string = '.'
        #     else:
        #         unit_string = 'followed by a space with a unit from the following list: {self.units}.'
            
        #     message = (f'The value, "{data}", for the subsection, "{subsection}", in the "{section}" section '
        #                f'should be a number{range_string}{unit_string}')
        #     return message
            
        #     return create_error_message(self.kwargs['section'], 
        #                                 self.kwargs['subsection'], 
        #                                 self.can_be_range, 
        #                                 self.kwargs['no_units'] if 'no_units' in self.kwargs else False)

    return NewCls



SubjectUnitChecker = partialclass(UnitChecker, error_message_builder = create_error_message2,  section = 'SUBJECT')
test = SubjectUnitChecker(['g'], True, subsection = 'SUBSECTION')
test.validate('5 - 6 g')


def test_message_builder(self, **kwargs):
    print(self)
    print(kwargs)

class TestClass(object):
    def __init__(self, message_builder, **kwargs):
        self.message_builder = message_builder
        self.kwargs = kwargs
    
    def _build_message(self):
        return self.message_builder(self)

test = TestClass(test_message_builder)
test._build_message()






metabolomics_workbench_schema = \
{'type': 'object',
 'properties': {'STUDY_ID': {'type': 'string'},
                'ANALYSIS_ID': {'type': 'string'},
                'VERSION': {'type': 'string'},
                'CREATED_ON': {'type': 'string'},
                'PROJECT_ID': {'type': 'string'},
                'HEADER': {'type': 'string'},
                'DATATRACK_ID': {'type': 'string'},
                'filename': {'type': 'string'}},
 'required': ['VERSION', 'CREATED_ON'],
 'additionalProperties': False}

project_schema = \
{'type': 'object',
 'properties': {'PROJECT_TITLE': {'type': 'string'},
                'PROJECT_TYPE': {'type': 'string'},
                'PROJECT_SUMMARY': {'type': 'string'},
                'INSTITUTE': {'type': 'string'},
                'DEPARTMENT': {'type': 'string'},
                'LABORATORY': {'type': 'string'},
                'LAST_NAME': {'type': 'string'},
                'FIRST_NAME': {'type': 'string'},
                'ADDRESS': {'type': 'string'},
                'EMAIL': {'type': 'string'},
                'PHONE': {'type': 'string'},
                'FUNDING_SOURCE': {'type': 'string'},
                'PROJECT_COMMENTS': {'type': 'string'},
                'PUBLICATIONS': {'type': 'string'},
                'CONTRIBUTORS': {'type': 'string'},
                'DOI': {'type': 'string'}},
 'required': ['PROJECT_TITLE',
              'PROJECT_SUMMARY',
              'INSTITUTE',
              'LAST_NAME',
              'FIRST_NAME',
              'ADDRESS',
              'EMAIL',
              'PHONE'],
 'additionalProperties': False}

study_schema = \
{'type': 'object',
 'properties': {'STUDY_TITLE': {'type': 'string'},
                'STUDY_TYPE': {'type': 'string'},
                'STUDY_SUMMARY': {'type': 'string'},
                'INSTITUTE': {'type': 'string'},
                'DEPARTMENT': {'type': 'string'},
                'LABORATORY': {'type': 'string'},
                'LAST_NAME': {'type': 'string'},
                'FIRST_NAME': {'type': 'string'},
                'ADDRESS': {'type': 'string'},
                'EMAIL': {'type': 'string'},
                'PHONE': {'type': 'string'},
                'SUBMIT_DATE': {'type': 'string'},
                'NUM_GROUPS': {'type': 'string'},
                'TOTAL_SUBJECTS': {'type': 'string'},
                'NUM_MALES': {'type': 'string'},
                'NUM_FEMALES': {'type': 'string'},
                'STUDY_COMMENTS': {'type': 'string'},
                'PUBLICATIONS': {'type': 'string'}},
 'required': ['STUDY_TITLE',
              'STUDY_SUMMARY',
              'INSTITUTE',
              'LAST_NAME',
              'FIRST_NAME',
              'ADDRESS',
              'EMAIL',
              'PHONE'],
 'additionalProperties': False}

subject_schema = \
{'type': 'object',
 'properties': {'SUBJECT_TYPE': {'type': 'string'},
                'SUBJECT_SPECIES': {'type': 'string'},
                'TAXONOMY_ID': {'type': 'string'},
                'GENOTYPE_STRAIN': {'type': 'string'},
                'AGE_OR_AGE_RANGE': {'type': 'string', 'pattern': create_units_regex(['weeks', 'days', 'months', 'years'], True)},
                'WEIGHT_OR_WEIGHT_RANGE': {'type': 'string', 'pattern': create_units_regex(['g', 'mg', 'kg', 'lbs'], True)},
                'HEIGHT_OR_HEIGHT_RANGE': {'type': 'string', 'pattern': create_units_regex(['cm', 'in'], True)},
                'GENDER': {'type': 'string'},
                'HUMAN_RACE': {'type': 'string'},
                'HUMAN_ETHNICITY': {'type': 'string'},
                'HUMAN_TRIAL_TYPE': {'type': 'string'},
                'HUMAN_LIFESTYLE_FACTORS': {'type': 'string'},
                'HUMAN_MEDICATIONS': {'type': 'string'},
                'HUMAN_PRESCRIPTION_OTC': {'type': 'string'},
                'HUMAN_SMOKING_STATUS': {'type': 'string'},
                'HUMAN_ALCOHOL_DRUG_USE': {'type': 'string'},
                'HUMAN_NUTRITION': {'type': 'string'},
                'HUMAN_INCLUSION_CRITERIA': {'type': 'string'},
                'HUMAN_EXCLUSION_CRITERIA': {'type': 'string'},
                'ANIMAL_ANIMAL_SUPPLIER': {'type': 'string'},
                'ANIMAL_HOUSING': {'type': 'string'},
                'ANIMAL_LIGHT_CYCLE': {'type': 'string'},
                'ANIMAL_FEED': {'type': 'string'},
                'ANIMAL_WATER': {'type': 'string'},
                'ANIMAL_INCLUSION_CRITERIA': {'type': 'string'},
                'CELL_BIOSOURCE_OR_SUPPLIER': {'type': 'string'},
                'CELL_STRAIN_DETAILS': {'type': 'string'},
                'SUBJECT_COMMENTS': {'type': 'string'},
                'CELL_PRIMARY_IMMORTALIZED': {'type': 'string'},
                'CELL_PASSAGE_NUMBER': {'type': 'string'},
                'CELL_COUNTS': {'type': 'string'},
                'SPECIES_GROUP': {'type': 'string'}},
 'required': ['SUBJECT_TYPE', 'SUBJECT_SPECIES'],
 'additionalProperties': False}

subject_sample_factors_schema = \
{'type': 'array',
 'items': {'type': 'object',
           'properties': {'Subject ID': {'type': 'string'},
                          'Sample ID': {'type': 'string'},
                          'Factors': {'type': 'object'},
                          'Additional sample data': {'type': 'object',
                                                     'properties': {'RAW_FILE_NAME': {'type': 'string'}},
                                                     'additionalProperties': True}},
           'required': ['Subject ID', 'Sample ID', 'Factors'],
           'additionalProperties': False}}    

collection_schema = \
{'type': 'object',
 'properties': {'COLLECTION_SUMMARY': {'type': 'string'},
                'COLLECTION_PROTOCOL_ID': {'type': 'string'},
                'COLLECTION_PROTOCOL_FILENAME': {'type': 'string'},
                'COLLECTION_PROTOCOL_COMMENTS': {'type': 'string'},
                'SAMPLE_TYPE': {'type': 'string'},
                'COLLECTION_METHOD': {'type': 'string'},
                'COLLECTION_LOCATION': {'type': 'string'},
                'COLLECTION_FREQUENCY': {'type': 'string'},
                'COLLECTION_DURATION': {'type': 'string'},
                'COLLECTION_TIME': {'type': 'string'},
                'VOLUMEORAMOUNT_COLLECTED': {'type': 'string'},
                'STORAGE_CONDITIONS': {'type': 'string'},
                'COLLECTION_VIALS': {'type': 'string'},
                'STORAGE_VIALS': {'type': 'string'},
                'COLLECTION_TUBE_TEMP': {'type': 'string'},
                'ADDITIVES': {'type': 'string'},
                'BLOOD_SERUM_OR_PLASMA': {'type': 'string', 'enum':['Blood', 'Serum', 'Plasma']},
                'TISSUE_CELL_IDENTIFICATION': {'type': 'string'},
                'TISSUE_CELL_QUANTITY_TAKEN': {'type': 'string'}},
 'required': ['COLLECTION_SUMMARY'],
 'additionalProperties': False}

treatment_schema = \
{'type': 'object',
 'properties': {'TREATMENT_SUMMARY': {'type': 'string'},
                'TREATMENT_PROTOCOL_ID': {'type': 'string'},
                'TREATMENT_PROTOCOL_FILENAME': {'type': 'string'},
                'TREATMENT_PROTOCOL_COMMENTS': {'type': 'string'},
                'TREATMENT': {'type': 'string'},
                'TREATMENT_COMPOUND': {'type': 'string'},
                'TREATMENT_ROUTE': {'type': 'string'},
                'TREATMENT_DOSE': {'type': 'string'},
                'TREATMENT_DOSEVOLUME': {'type': 'string'},
                'TREATMENT_DOSEDURATION': {'type': 'string', 'pattern': create_units_regex(['h', 'weeks', 'days'], True)},
                'TREATMENT_VEHICLE': {'type': 'string'},
                'ANIMAL_VET_TREATMENTS': {'type': 'string'},
                'ANIMAL_ANESTHESIA': {'type': 'string'},
                'ANIMAL_ACCLIMATION_DURATION': {'type': 'string'},
                'ANIMAL_FASTING': {'type': 'string'},
                'ANIMAL_ENDP_EUTHANASIA': {'type': 'string'},
                'ANIMAL_ENDP_TISSUE_COLL_LIST': {'type': 'string'},
                'ANIMAL_ENDP_TISSUE_PROC_METHOD': {'type': 'string'},
                'ANIMAL_ENDP_CLINICAL_SIGNS': {'type': 'string'},
                'HUMAN_FASTING': {'type': 'string'},
                'HUMAN_ENDP_CLINICAL_SIGNS': {'type': 'string'},
                'CELL_STORAGE': {'type': 'string'},
                'CELL_GROWTH_CONTAINER': {'type': 'string'},
                'CELL_GROWTH_CONFIG': {'type': 'string'},
                'CELL_GROWTH_RATE': {'type': 'string'},
                'CELL_INOC_PROC': {'type': 'string'},
                'CELL_MEDIA': {'type': 'string'},
                'CELL_ENVIR_COND': {'type': 'string'},
                'CELL_HARVESTING': {'type': 'string'},
                'PLANT_GROWTH_SUPPORT': {'type': 'string'},
                'PLANT_GROWTH_LOCATION': {'type': 'string'},
                'PLANT_PLOT_DESIGN': {'type': 'string'},
                'PLANT_LIGHT_PERIOD': {'type': 'string'},
                'PLANT_HUMIDITY': {'type': 'string'},
                'PLANT_TEMP': {'type': 'string', 'pattern': create_units_regex(['°C', 'C'], True)},
                'PLANT_WATERING_REGIME': {'type': 'string'},
                'PLANT_NUTRITIONAL_REGIME': {'type': 'string'},
                'PLANT_ESTAB_DATE': {'type': 'string'},
                'PLANT_HARVEST_DATE': {'type': 'string'},
                'PLANT_GROWTH_STAGE': {'type': 'string'},
                'PLANT_METAB_QUENCH_METHOD': {'type': 'string'},
                'PLANT_HARVEST_METHOD': {'type': 'string'},
                'PLANT_STORAGE': {'type': 'string'},
                'CELL_PCT_CONFLUENCE': {'type': 'string'},
                'CELL_MEDIA_LASTCHANGED': {'type': 'string'}},
 'required': ['TREATMENT_SUMMARY'],
 'additionalProperties': False}

sampleprep_schema = \
{'type': 'object',
 'properties': {'SAMPLEPREP_SUMMARY': {'type': 'string'},
                'SAMPLEPREP_PROTOCOL_ID': {'type': 'string'},
                'SAMPLEPREP_PROTOCOL_FILENAME': {'type': 'string'},
                'SAMPLEPREP_PROTOCOL_COMMENTS': {'type': 'string'},
                'PROCESSING_METHOD': {'type': 'string'},
                'PROCESSING_STORAGE_CONDITIONS': {'type': 'string'},
                'EXTRACTION_METHOD': {'type': 'string'},
                'EXTRACT_CONCENTRATION_DILUTION': {'type': 'string'},
                'EXTRACT_ENRICHMENT': {'type': 'string'},
                'EXTRACT_CLEANUP': {'type': 'string'},
                'EXTRACT_STORAGE': {'type': 'string'},
                'SAMPLE_RESUSPENSION': {'type': 'string'},
                'SAMPLE_DERIVATIZATION': {'type': 'string'},
                'SAMPLE_SPIKING': {'type': 'string'},
                'ORGAN': {'type': 'string'},
                'ORGAN_SPECIFICATION': {'type': 'string'},
                'CELL_TYPE': {'type': 'string'},
                'SUBCELLULAR_LOCATION': {'type': 'string'}},
 'required': ['SAMPLEPREP_SUMMARY'],
 'additionalProperties': False}

chromatography_schema = \
{'type': 'object',
 'properties': {'CHROMATOGRAPHY_SUMMARY': {'type': 'string'},
                'CHROMATOGRAPHY_TYPE': {'type': 'string'},
                'INSTRUMENT_NAME': {'type': 'string'},
                'COLUMN_NAME': {'type': 'string'},
                'FLOW_GRADIENT': {'type': 'string'},
                'FLOW_RATE': {'type': 'string', 'pattern': create_units_regex(['mL/min', 'uL/min', 'μL/min'], True)},
                'COLUMN_TEMPERATURE': {'type': 'string', 'pattern': create_units_regex(['°C', 'C'], True)},
                'METHODS_FILENAME': {'type': 'string'},
                'SAMPLE_INJECTION': {'type': 'string', 'pattern': create_units_regex(['μL', 'uL'])},
                'SOLVENT_A': {'type': 'string'},
                'SOLVENT_B': {'type': 'string'},
                'METHODS_ID': {'type': 'string'},
                'COLUMN_PRESSURE': {'type': 'string', 'pattern': create_units_regex(['psi', 'bar'], True)},
                'INJECTION_TEMPERATURE': {'type': 'string', 'pattern': create_units_regex(['°C', 'C'], True)},
                'INTERNAL_STANDARD': {'type': 'string'},
                'INTERNAL_STANDARD_MT': {'type': 'string'},
                'RETENTION_INDEX': {'type': 'string'},
                'RETENTION_TIME': {'type': 'string'},
                'SAMPLING_CONE': {'type': 'string'},
                'ANALYTICAL_TIME': {'type': 'string', 'pattern': create_units_regex(['min'], True)},
                'CAPILLARY_VOLTAGE': {'type': 'string', 'pattern': create_units_regex(['V', 'kV'])},
                'MIGRATION_TIME': {'type': 'string'},
                'OVEN_TEMPERATURE': {'type': 'string'},
                'PRECONDITIONING': {'type': 'string'},
                'RUNNING_BUFFER': {'type': 'string'},
                'RUNNING_VOLTAGE': {'type': 'string', 'pattern': create_units_regex(['V', 'kV'])},
                'SHEATH_LIQUID': {'type': 'string'},
                'TIME_PROGRAM': {'type': 'string'},
                'TRANSFERLINE_TEMPERATURE': {'type': 'string', 'pattern': create_units_regex(['°C', 'C'])},
                'WASHING_BUFFER': {'type': 'string'},
                'WEAK_WASH_SOLVENT_NAME': {'type': 'string'},
                'WEAK_WASH_VOLUME': {'type': 'string', 'pattern': create_units_regex(['μL', 'uL'])},
                'STRONG_WASH_SOLVENT_NAME': {'type': 'string'},
                'STRONG_WASH_VOLUME': {'type': 'string', 'pattern': create_units_regex(['μL', 'uL'])},
                'TARGET_SAMPLE_TEMPERATURE': {'type': 'string', 'pattern': create_units_regex(['°C', 'C'])},
                'SAMPLE_LOOP_SIZE': {'type': 'string', 'pattern': create_units_regex(['μL', 'uL'])},
                'SAMPLE_SYRINGE_SIZE': {'type': 'string', 'pattern': create_units_regex(['μL', 'uL'])},
                'RANDOMIZATION_ORDER': {'type': 'string'},
                'CHROMATOGRAPHY_COMMENTS': {'type': 'string'}},
 'required': ['CHROMATOGRAPHY_TYPE',
              'INSTRUMENT_NAME',
              'COLUMN_NAME',
              'FLOW_GRADIENT',
              'FLOW_RATE',
              'COLUMN_TEMPERATURE',
              'SOLVENT_A',
              'SOLVENT_B'],
 'additionalProperties': False}

analysis_schema = \
{'type': 'object',
 'properties': {'ANALYSIS_TYPE': {'type': 'string', 'enum': ['MS', 'NMR']},
                'LABORATORY_NAME': {'type': 'string'},
                'ACQUISITION_DATE': {'type': 'string'},
                'SOFTWARE_VERSION': {'type': 'string'},
                'OPERATOR_NAME': {'type': 'string'},
                'DETECTOR_TYPE': {'type': 'string'},
                'ANALYSIS_PROTOCOL_FILE': {'type': 'string'},
                'ACQUISITION_PARAMETERS_FILE': {'type': 'string'},
                'PROCESSING_PARAMETERS_FILE': {'type': 'string'},
                'DATA_FORMAT': {'type': 'string'},
                'ACQUISITION_ID': {'type': 'string'},
                'ACQUISITION_TIME': {'type': 'string'},
                'ANALYSIS_COMMENTS': {'type': 'string'},
                'ANALYSIS_DISPLAY': {'type': 'string'},
                'INSTRUMENT_NAME': {'type': 'string'},
                'INSTRUMENT_PARAMETERS_FILE': {'type': 'string'},
                'NUM_FACTORS': {'type': 'string'},
                'NUM_METABOLITES': {'type': 'string'},
                'PROCESSED_FILE': {'type': 'string'},
                'RANDOMIZATION_ORDER': {'type': 'string'},
                'RAW_FILE': {'type': 'string'}},
 'required': ['ANALYSIS_TYPE'],
 'additionalProperties': False}

results_file_schema = \
{'type': 'object',
 'properties': {'filename': {'type': 'string'},
  'UNITS': {'type': 'string'},
  'Has m/z': {'type': 'string'},
  'Has RT': {'type': 'string'},
  'RT units': {'type': 'string'}},
 'required': ['filename', 'UNITS'],
 'additionalProperties': False}

ms_schema = \
{'type': 'object',
 'properties': {'INSTRUMENT_NAME': {'type': 'string'},
                'INSTRUMENT_TYPE': {'type': 'string'},
                'MS_TYPE': {'type': 'string'},
                'ION_MODE': {'type': 'string'},
                'CAPILLARY_TEMPERATURE': {'type': 'string', 'pattern': create_units_regex(['°C', 'C'], True)},
                'CAPILLARY_VOLTAGE': {'type': 'string', 'pattern': create_units_regex(['V', 'kV'])},
                'COLLISION_ENERGY': {'type': 'string'},
                'COLLISION_GAS': {'type': 'string', 'enum': ['Nitrogen', 'Argon']},
                'DRY_GAS_FLOW': {'type': 'string', 'pattern': create_units_regex(['L/hr', 'L/min'])},
                'DRY_GAS_TEMP': {'type': 'string', 'pattern': create_units_regex(['°C', 'C'])},
                'FRAGMENT_VOLTAGE': {'type': 'string', 'pattern': create_units_regex(['V'])},
                'FRAGMENTATION_METHOD': {'type': 'string'},
                'GAS_PRESSURE': {'type': 'string', 'pattern': create_units_regex(['psi', 'psig', 'bar', 'kPa'])},
                'HELIUM_FLOW': {'type': 'string', 'pattern': create_units_regex(['mL/min'])},
                'ION_SOURCE_TEMPERATURE': {'type': 'string', 'pattern': create_units_regex(['°C', 'C'])},
                'ION_SPRAY_VOLTAGE': {'type': 'string', 'pattern': create_units_regex(['V', 'kV'])},
                'IONIZATION': {'type': 'string'},
                'IONIZATION_ENERGY': {'type': 'string', 'pattern': create_units_regex(['eV'])},
                'IONIZATION_POTENTIAL': {'type': 'string'},
                'MASS_ACCURACY': {'type': 'string'},
                'PRECURSOR_TYPE': {'type': 'string'},
                'REAGENT_GAS': {'type': 'string'},
                'SOURCE_TEMPERATURE': {'type': 'string', 'pattern': create_units_regex(['°C', 'C'])},
                'SPRAY_VOLTAGE': {'type': 'string', 'pattern': create_units_regex(['kV'])},
                'ACTIVATION_PARAMETER': {'type': 'string'},
                'ACTIVATION_TIME': {'type': 'string', 'pattern': create_units_regex(['ms'])},
                'ATOM_GUN_CURRENT': {'type': 'string'},
                'AUTOMATIC_GAIN_CONTROL': {'type': 'string'},
                'BOMBARDMENT': {'type': 'string'},
                'CDL_SIDE_OCTOPOLES_BIAS_VOLTAGE': {'type': 'string'},
                'CDL_TEMPERATURE': {'type': 'string'},
                'DATAFORMAT': {'type': 'string'},
                'DESOLVATION_GAS_FLOW': {'type': 'string', 'pattern': create_units_regex(['L/hr', 'L/min'])},
                'DESOLVATION_TEMPERATURE': {'type': 'string', 'pattern': create_units_regex(['°C', 'C'])},
                'INTERFACE_VOLTAGE': {'type': 'string'},
                'IT_SIDE_OCTOPOLES_BIAS_VOLTAGE': {'type': 'string'},
                'LASER': {'type': 'string'},
                'MATRIX': {'type': 'string'},
                'NEBULIZER': {'type': 'string'},
                'OCTPOLE_VOLTAGE': {'type': 'string', 'pattern': create_units_regex(['V'])},
                'PROBE_TIP': {'type': 'string'},
                'RESOLUTION_SETTING': {'type': 'string'},
                'SAMPLE_DRIPPING': {'type': 'string'},
                'SCAN_RANGE_MOVERZ': {'type': 'string'},
                'SCANNING': {'type': 'string'},
                'SCANNING_CYCLE': {'type': 'string'},
                'SCANNING_RANGE': {'type': 'string'},
                'SKIMMER_VOLTAGE': {'type': 'string', 'pattern': create_units_regex(['V'])},
                'TUBE_LENS_VOLTAGE': {'type': 'string'},
                'MS_COMMENTS': {'type': 'string'},
                'MS_RESULTS_FILE': results_file_schema},
 'required': ['INSTRUMENT_NAME', 'INSTRUMENT_TYPE', 'MS_TYPE', 'ION_MODE'],
 'additionalProperties': False}

nmr_schema = \
{'type': 'object',
 'properties': {'INSTRUMENT_NAME': {'type': 'string'},
                'INSTRUMENT_TYPE': {'type': 'string'},
                'NMR_EXPERIMENT_TYPE': {'type': 'string'},
                'NMR_COMMENTS': {'type': 'string'},
                'FIELD_FREQUENCY_LOCK': {'type': 'string'},
                'STANDARD_CONCENTRATION': {'type': 'string', 'pattern': create_units_regex(['mM'])},
                'SPECTROMETER_FREQUENCY': {'type': 'string', 'pattern': create_units_regex(['MHz'])},
                'NMR_PROBE': {'type': 'string'},
                'NMR_SOLVENT': {'type': 'string'},
                'NMR_TUBE_SIZE': {'type': 'string'},
                'SHIMMING_METHOD': {'type': 'string'},
                'PULSE_SEQUENCE': {'type': 'string'},
                'WATER_SUPPRESSION': {'type': 'string'},
                'PULSE_WIDTH': {'type': 'string'},
                'POWER_LEVEL': {'type': 'string', 'pattern': create_units_regex(['W'])},
                'RECEIVER_GAIN': {'type': 'string', 'pattern': '^\d+$'},
                'OFFSET_FREQUENCY': {'type': 'string', 'pattern': create_units_regex(['ppm', 'Hz'])},
                'PRESATURATION_POWER_LEVEL': {'type': 'string', 'pattern': create_units_regex(['W'])},
                'CHEMICAL_SHIFT_REF_CPD': {'type': 'string'},
                'TEMPERATURE': {'type': 'string', 'pattern': create_units_regex(['°C', 'C', 'K'])},
                'NUMBER_OF_SCANS': {'type': 'string', 'pattern': '^\d+$'},
                'DUMMY_SCANS': {'type': 'string'},
                'ACQUISITION_TIME': {'type': 'string', 'pattern': create_units_regex(['s'])},
                'RELAXATION_DELAY': {'type': 'string', 'pattern': create_units_regex(['s', 'ms', 'us', 'μs'])},
                'SPECTRAL_WIDTH': {'type': 'string', 'pattern': create_units_regex(['ppm', 'Hz'])},
                'NUM_DATA_POINTS_ACQUIRED': {'type': 'string', 'pattern': '^\d+$'},
                'REAL_DATA_POINTS': {'type': 'string'},
                'LINE_BROADENING': {'type': 'string', 'pattern': create_units_regex(['Hz'])},
                'ZERO_FILLING': {'type': 'string'},
                'APODIZATION': {'type': 'string'},
                'BASELINE_CORRECTION_METHOD': {'type': 'string'},
                'CHEMICAL_SHIFT_REF_STD': {'type': 'string'},
                'BINNED_INCREMENT': {'type': 'string', 'pattern': create_units_regex(['ppm'])},
                'BINNED_DATA_NORMALIZATION_METHOD': {'type': 'string'},
                'BINNED_DATA_PROTOCOL_FILE': {'type': 'string'},
                'BINNED_DATA_CHEMICAL_SHIFT_RANGE': {'type': 'string'},
                'BINNED_DATA_EXCLUDED_RANGE': {'type': 'string'},
                'NMR_RESULTS_FILE': results_file_schema},
 'required': ['INSTRUMENT_NAME',
              'INSTRUMENT_TYPE',
              'NMR_EXPERIMENT_TYPE',
              'SPECTROMETER_FREQUENCY'],
 'additionalProperties': False}

data_schema = \
{'type': 'array',
 'items': {'type': 'object',
           'properties': {'Metabolite': {'type': 'string'},
                          'Bin range(ppm)': {'type': 'string'}},
           'required': [],
           'additionalProperties': True}}

extended_schema = \
{'type': 'array',
 'items': {'type': 'object',
           'properties': {'Metabolite': {'type': 'string'},
                          'sample_id': {'type': 'string'}},
           'required': ['Metabolite', 'sample_id'],
           'additionalProperties': True}}

ms_metabolite_data_schema = \
{'type': 'object',
 'properties': {'Units': {'type': 'string'},
                'Data': data_schema,
                'Metabolites': data_schema,
                'Extended': extended_schema},
 'required': ['Units', 'Data'],
 'additionalProperties': False}

nmr_binned_data_schema = \
{'type': 'object',
 'properties': {'Units': {'type': 'string'},
                'Data': data_schema},
 'required': ['Units', 'Data'],
 'additionalProperties': False}

base_required_schema = \
{'properties': {'METABOLOMICS WORKBENCH': metabolomics_workbench_schema,
               'PROJECT': project_schema,
               'STUDY': study_schema,
               'SUBJECT': subject_schema,
               'SUBJECT_SAMPLE_FACTORS': subject_sample_factors_schema,
               'COLLECTION': collection_schema,
               'TREATMENT': treatment_schema,
               'SAMPLEPREP': sampleprep_schema,
               'ANALYSIS': analysis_schema},
 'required': ['METABOLOMICS WORKBENCH',
              'PROJECT',
              'STUDY',
              'SUBJECT',
              'SUBJECT_SAMPLE_FACTORS',
              'COLLECTION',
              'TREATMENT',
              'SAMPLEPREP',
              'ANALYSIS'],
 'additionalProperties': False}

ms_required_schema = deepcopy(base_required_schema)
ms_required_schema['properties']['MS'] = ms_schema
ms_required_schema['properties']['MS_METABOLITE_DATA'] = ms_metabolite_data_schema
ms_required_schema['properties']['CHROMATOGRAPHY'] = chromatography_schema
ms_required_schema['required'].extend(['MS'])
ms_required_schema['if'] = {'properties': {'MS_METABOLITE_DATA':{'not':{}}}}
ms_required_schema['then'] = {'properties': {'MS':{'required':['MS_RESULTS_FILE']}}}
# TODO catch the "'MS_RESULTS_FILE' is a required property" error and say that either a data section or results file is required.

nmr_required_schema = deepcopy(base_required_schema)
nmr_required_schema['properties']['NM'] = nmr_schema
nmr_required_schema['properties']['NMR_METABOLITE_DATA'] = ms_metabolite_data_schema
nmr_required_schema['properties']['NMR_BINNED_DATA'] = nmr_binned_data_schema
nmr_required_schema['required'].extend(['NM'])
nmr_required_schema['if'] = {'allOf':[{'properties': {'NMR_METABOLITE_DATA':{'not':{}}}}, {'properties': {'NMR_BINNED_DATA':{'not':{}}}}]}
nmr_required_schema['then'] = {'properties': {'NM':{'required':['NMR_RESULTS_FILE']}}}
# TODO catch the "'NMR_RESULTS_FILE' is a required property" error and say that either a data section or results file is required.


compiled_schema = \
{'type': 'object',
 'oneOf':[ms_required_schema,
          nmr_required_schema]}
    
    
# {'type': 'object',
#  'properties': {'METABOLOMICS WORKBENCH': metabolomics_workbench_schema,
#                 'PROJECT': project_schema,
#                 'STUDY': study_schema,
#                 'SUBJECT': subject_schema,
#                 'SUBJECT_SAMPLE_FACTORS': subject_sample_factors_schema,
#                 'COLLECTION': collection_schema,
#                 'TREATMENT': treatment_schema,
#                 'SAMPLEPREP': sampleprep_schema,
#                 'ANALYSIS': analysis_schema,
#                 'CHROMATOGRAPHY': chromatography_schema,
#                 'MS_METABOLITE_DATA': ms_metabolite_data_schema,
#                 'NMR_METABOLITE_DATA': ms_metabolite_data_schema,
#                 'NMR_BINNED_DATA': nmr_binned_data_schema,
#                 'MS': ms_schema,
#                 'NM': nmr_schema},
#  'required': ['METABOLOMICS WORKBENCH',
#               'PROJECT',
#               'STUDY',
#               'SUBJECT',
#               'SUBJECT_SAMPLE_FACTORS',
#               'COLLECTION',
#               'TREATMENT',
#               'SAMPLEPREP',
#               'ANALYSIS'],
#  'dependentRequired': {},
#  'additionalProperties': False}










metabolomics_workbench_schema = Schema(
    {
        Optional("STUDY_ID"): str,
        Optional("ANALYSIS_ID"): str,
        "VERSION": str,
        "CREATED_ON": str,
        Optional("PROJECT_ID"): str,
        Optional("HEADER"): str,
        Optional("DATATRACK_ID"): str,
        Optional("filename"): str
    }
)

project_schema = Schema(
    {
        "PROJECT_TITLE": str,
        Optional("PROJECT_TYPE"): str,
        "PROJECT_SUMMARY": str,
        "INSTITUTE": str,
        Optional("DEPARTMENT"): str,
        Optional("LABORATORY"): str,
        "LAST_NAME": str,
        "FIRST_NAME": str,
        "ADDRESS": str,
        "EMAIL": str,
        "PHONE": str,
        Optional("FUNDING_SOURCE"): str,
        Optional("PROJECT_COMMENTS"): str,
        Optional("PUBLICATIONS"): str,
        Optional("CONTRIBUTORS"): str,
        Optional("DOI"): str
    }
)

study_schema = Schema(
    {
        "STUDY_TITLE": str,
        Optional("STUDY_TYPE"): str,
        "STUDY_SUMMARY": str,
        "INSTITUTE": str,
        Optional("DEPARTMENT"): str,
        Optional("LABORATORY"): str,
        "LAST_NAME": str,
        "FIRST_NAME": str,
        "ADDRESS": str,
        "EMAIL": str,
        "PHONE": str,
        Optional("SUBMIT_DATE"): str,  # assumed
        Optional("NUM_GROUPS"): str,
        Optional("TOTAL_SUBJECTS"): str,
        Optional("NUM_MALES"): str,
        Optional("NUM_FEMALES"): str,
        Optional("STUDY_COMMENTS"): str,
        Optional("PUBLICATIONS"): str
    }
)

subject_error_message = partial(create_error_message, 'SUBJECT')
subject_schema = Schema(
    {
        "SUBJECT_TYPE": str,
        "SUBJECT_SPECIES": str,
        Optional("TAXONOMY_ID"): str,
        Optional("GENOTYPE_STRAIN"): str,
        Optional("AGE_OR_AGE_RANGE"): And(str, UnitChecker(['weeks', 'days', 'months', 'years'], True, subject_error_message('AGE_OR_AGE_RANGE', True))),
        Optional("WEIGHT_OR_WEIGHT_RANGE"): And(str, UnitChecker(['g', 'mg', 'kg', 'lbs'], True, subject_error_message('WEIGHT_OR_WEIGHT_RANGE', True))),
        Optional("HEIGHT_OR_HEIGHT_RANGE"): And(str, UnitChecker(['cm', 'in'], True, subject_error_message('HEIGHT_OR_HEIGHT_RANGE', True))),
        Optional("GENDER"): str,
        Optional("HUMAN_RACE"): str,
        Optional("HUMAN_ETHNICITY"): str,
        Optional("HUMAN_TRIAL_TYPE"): str,
        Optional("HUMAN_LIFESTYLE_FACTORS"): str,
        Optional("HUMAN_MEDICATIONS"): str,
        Optional("HUMAN_PRESCRIPTION_OTC"): str,
        Optional("HUMAN_SMOKING_STATUS"): str,
        Optional("HUMAN_ALCOHOL_DRUG_USE"): str,
        Optional("HUMAN_NUTRITION"): str,
        Optional("HUMAN_INCLUSION_CRITERIA"): str,
        Optional("HUMAN_EXCLUSION_CRITERIA"): str,
        Optional("ANIMAL_ANIMAL_SUPPLIER"): str,
        Optional("ANIMAL_HOUSING"): str,
        Optional("ANIMAL_LIGHT_CYCLE"): str,
        Optional("ANIMAL_FEED"): str,
        Optional("ANIMAL_WATER"): str,
        Optional("ANIMAL_INCLUSION_CRITERIA"): str,
        Optional("CELL_BIOSOURCE_OR_SUPPLIER"): str,
        Optional("CELL_STRAIN_DETAILS"): str,
        Optional("SUBJECT_COMMENTS"): str,
        Optional("CELL_PRIMARY_IMMORTALIZED"): str,
        Optional("CELL_PASSAGE_NUMBER"): str,
        Optional("CELL_COUNTS"): str,
        Optional("SPECIES_GROUP"): str
    }
)

subject_sample_factors_schema = Schema(
    [
        {
            "Subject ID": str,
            "Sample ID": str,
            "Factors": dict,
            Optional("Additional sample data"): {
                Optional("RAW_FILE_NAME"): str,
                Optional(str): str
            }
        }
    ]
)

collection_schema = Schema(
    {
        "COLLECTION_SUMMARY": str,
        Optional("COLLECTION_PROTOCOL_ID"): str,
        Optional("COLLECTION_PROTOCOL_FILENAME"): str,
        Optional("COLLECTION_PROTOCOL_COMMENTS"): str,
        Optional("SAMPLE_TYPE"): str,  # required as of mwTab file format specification 1.5
        Optional("COLLECTION_METHOD"): str,
        Optional("COLLECTION_LOCATION"): str,
        Optional("COLLECTION_FREQUENCY"): str,
        Optional("COLLECTION_DURATION"): str,
        Optional("COLLECTION_TIME"): str,
        Optional("VOLUMEORAMOUNT_COLLECTED"): str,
        Optional("STORAGE_CONDITIONS"): str,
        Optional("COLLECTION_VIALS"): str,
        Optional("STORAGE_VIALS"): str,
        Optional("COLLECTION_TUBE_TEMP"): str,
        Optional("ADDITIVES"): str,
        # TODO Ask Hunter if these should be lowercase or leave them as is.
        Optional("BLOOD_SERUM_OR_PLASMA"): And(str, SetChecker(['Blood', 'Serum', 'Plasma'])),
        Optional("TISSUE_CELL_IDENTIFICATION"): str,
        Optional("TISSUE_CELL_QUANTITY_TAKEN"): str
    }
)

treatment_error_message = partial(create_error_message, 'TREATMENT')
treatment_schema = Schema(
    {
        "TREATMENT_SUMMARY": str,
        Optional("TREATMENT_PROTOCOL_ID"): str,
        Optional("TREATMENT_PROTOCOL_FILENAME"): str,
        Optional("TREATMENT_PROTOCOL_COMMENTS"): str,
        Optional("TREATMENT"): str,
        Optional("TREATMENT_COMPOUND"): str,
        Optional("TREATMENT_ROUTE"): str,
        Optional("TREATMENT_DOSE"): str,
        Optional("TREATMENT_DOSEVOLUME"): str,
        Optional("TREATMENT_DOSEDURATION"): And(str, UnitChecker(['h', 'days', 'weeks'], True, treatment_error_message('TREATMENT_DOSEDURATION', True))),
        Optional("TREATMENT_VEHICLE"): str,
        Optional("ANIMAL_VET_TREATMENTS"): str,
        Optional("ANIMAL_ANESTHESIA"): str,
        Optional("ANIMAL_ACCLIMATION_DURATION"): str,
        Optional("ANIMAL_FASTING"): str,
        Optional("ANIMAL_ENDP_EUTHANASIA"): str,
        Optional("ANIMAL_ENDP_TISSUE_COLL_LIST"): str,
        Optional("ANIMAL_ENDP_TISSUE_PROC_METHOD"): str,
        Optional("ANIMAL_ENDP_CLINICAL_SIGNS"): str,
        Optional("HUMAN_FASTING"): str,
        Optional("HUMAN_ENDP_CLINICAL_SIGNS"): str,
        Optional("CELL_STORAGE"): str,
        Optional("CELL_GROWTH_CONTAINER"): str,
        Optional("CELL_GROWTH_CONFIG"): str,
        Optional("CELL_GROWTH_RATE"): str,
        Optional("CELL_INOC_PROC"): str,
        Optional("CELL_MEDIA"): str,
        Optional("CELL_ENVIR_COND"): str,
        Optional("CELL_HARVESTING"): str,
        Optional("PLANT_GROWTH_SUPPORT"): str,
        Optional("PLANT_GROWTH_LOCATION"): str,
        Optional("PLANT_PLOT_DESIGN"): str,
        Optional("PLANT_LIGHT_PERIOD"): str,
        Optional("PLANT_HUMIDITY"): str,
        Optional("PLANT_TEMP"): And(str, UnitChecker(['°C', 'C'], True, treatment_error_message('PLANT_TEMP', True))),
        Optional("PLANT_WATERING_REGIME"): str,
        Optional("PLANT_NUTRITIONAL_REGIME"): str,
        Optional("PLANT_ESTAB_DATE"): str,
        Optional("PLANT_HARVEST_DATE"): str,
        Optional("PLANT_GROWTH_STAGE"): str,
        Optional("PLANT_METAB_QUENCH_METHOD"): str,
        Optional("PLANT_HARVEST_METHOD"): str,
        Optional("PLANT_STORAGE"): str,
        Optional("CELL_PCT_CONFLUENCE"): str,
        Optional("CELL_MEDIA_LASTCHANGED"): str
    }
)

sampleprep_schema = Schema(
    {
        "SAMPLEPREP_SUMMARY": str,
        Optional("SAMPLEPREP_PROTOCOL_ID"): str,
        Optional("SAMPLEPREP_PROTOCOL_FILENAME"): str,
        Optional("SAMPLEPREP_PROTOCOL_COMMENTS"): str,
        Optional("PROCESSING_METHOD"): str,
        Optional("PROCESSING_STORAGE_CONDITIONS"): str,
        Optional("EXTRACTION_METHOD"): str,
        Optional("EXTRACT_CONCENTRATION_DILUTION"): str,
        Optional("EXTRACT_ENRICHMENT"): str,
        Optional("EXTRACT_CLEANUP"): str,
        Optional("EXTRACT_STORAGE"): str,
        Optional("SAMPLE_RESUSPENSION"): str,
        Optional("SAMPLE_DERIVATIZATION"): str,
        Optional("SAMPLE_SPIKING"): str,
        Optional("ORGAN"): str,
        Optional("ORGAN_SPECIFICATION"): str,
        Optional("CELL_TYPE"): str,
        Optional("SUBCELLULAR_LOCATION"): str
    }
)

chromatography_error_message = partial(create_error_message, 'CHROMATOGRAPHY')
chromatography_schema = Schema(
    {
        Optional("CHROMATOGRAPHY_SUMMARY"): str,
        "CHROMATOGRAPHY_TYPE": str,
        "INSTRUMENT_NAME": str,
        "COLUMN_NAME": str,
        "FLOW_GRADIENT": str,
        # TODO ask Hunter about μ vs u. Allow both or force one?
        "FLOW_RATE": And(str, UnitChecker(['mL/min', 'uL/min', 'μL/min'], True, chromatography_error_message('FLOW_RATE', True))),
        # TODO ask Hunter about °C and C.
        "COLUMN_TEMPERATURE": And(str, UnitChecker(['°C', 'C'], True, chromatography_error_message('COLUMN_TEMPERATURE', True))),
        Optional("METHODS_FILENAME"): str,
        Optional("SAMPLE_INJECTION"): And(str, UnitChecker(['μL', 'uL'], error_message = chromatography_error_message('SAMPLE_INJECTION'))),
        "SOLVENT_A": str,
        "SOLVENT_B": str,
        Optional("METHODS_ID"): str,
        Optional("COLUMN_PRESSURE"): And(str, UnitChecker(['psi', 'bar'], True, chromatography_error_message('COLUMN_PRESSURE', True))),
        # TODO ask Hunter about a special case for "room temperature".
        Optional("INJECTION_TEMPERATURE"): And(str, UnitChecker(['°C', 'C'], True, chromatography_error_message('INJECTION_TEMPERATURE', True))),
        Optional("INTERNAL_STANDARD"): str,
        Optional("INTERNAL_STANDARD_MT"): str,
        Optional("RETENTION_INDEX"): str,
        Optional("RETENTION_TIME"): str,
        Optional("SAMPLING_CONE"): str,
        Optional("ANALYTICAL_TIME"): And(str, UnitChecker(['min'], True, chromatography_error_message('ANALYTICAL_TIME', True))),
        Optional("CAPILLARY_VOLTAGE"): And(str, UnitChecker(['V', 'kV'], error_message = chromatography_error_message('CAPILLARY_VOLTAGE'))),
        Optional("MIGRATION_TIME"): str,
        # TODO Ask hunter about cases like this one where something like "50°C for 1 min, then ramped at 20°C/min to 330°C, held constant for 5 min" is more common than a single number.
        Optional("OVEN_TEMPERATURE"): str,
        Optional("PRECONDITIONING"): str,
        Optional("RUNNING_BUFFER"): str,
        Optional("RUNNING_VOLTAGE"): And(str, UnitChecker(['V', 'kV'], error_message = chromatography_error_message('RUNNING_VOLTAGE'))),
        Optional("SHEATH_LIQUID"): str,
        Optional("TIME_PROGRAM"): str,
        Optional("TRANSFERLINE_TEMPERATURE"): And(str, UnitChecker(['°C', 'C'], error_message = chromatography_error_message('TRANSFERLINE_TEMPERATURE'))),
        Optional("WASHING_BUFFER"): str,
        Optional("WEAK_WASH_SOLVENT_NAME"): str,
        Optional("WEAK_WASH_VOLUME"): And(str, UnitChecker(['μL', 'uL'], error_message = chromatography_error_message('WEAK_WASH_VOLUME'))),
        Optional("STRONG_WASH_SOLVENT_NAME"): str,
        Optional("STRONG_WASH_VOLUME"): And(str, UnitChecker(['μL', 'uL'], error_message = chromatography_error_message('STRONG_WASH_VOLUME'))),
        Optional("TARGET_SAMPLE_TEMPERATURE"): And(str, UnitChecker(['°C', 'C'], error_message = chromatography_error_message('TARGET_SAMPLE_TEMPERATURE'))),
        Optional("SAMPLE_LOOP_SIZE"): And(str, UnitChecker(['μL', 'uL'], error_message = chromatography_error_message('SAMPLE_LOOP_SIZE'))),
        Optional("SAMPLE_SYRINGE_SIZE"): And(str, UnitChecker(['μL', 'uL'], error_message = chromatography_error_message('SAMPLE_SYRINGE_SIZE'))),
        Optional("RANDOMIZATION_ORDER"): str,
        Optional("CHROMATOGRAPHY_COMMENTS"): str
    }
)

analysis_schema = Schema(
    {
        "ANALYSIS_TYPE": And(str, SetChecker(['MS', 'NMR'])),
        Optional("LABORATORY_NAME"): str,
        Optional("ACQUISITION_DATE"): str,
        Optional("SOFTWARE_VERSION"): str,
        Optional("OPERATOR_NAME"): str,
        Optional("DETECTOR_TYPE"): str,
        Optional("ANALYSIS_PROTOCOL_FILE"): str,
        Optional("ACQUISITION_PARAMETERS_FILE"): str,
        Optional("PROCESSING_PARAMETERS_FILE"): str,
        Optional("DATA_FORMAT"): str,

        # not specified in mwTab specification (assumed)
        Optional("ACQUISITION_ID"): str,
        Optional("ACQUISITION_TIME"): str,
        Optional("ANALYSIS_COMMENTS"): str,
        Optional("ANALYSIS_DISPLAY"): str,
        Optional("INSTRUMENT_NAME"): str,
        Optional("INSTRUMENT_PARAMETERS_FILE"): str,
        Optional("NUM_FACTORS"): str,
        Optional("NUM_METABOLITES"): str,
        Optional("PROCESSED_FILE"): str,
        Optional("RANDOMIZATION_ORDER"): str,
        Optional("RAW_FILE"): str,
    }
)

results_file_schema = Schema(
    {
     'filename': str,
     'UNITS': str,
     Optional('Has m/z'): str,
     Optional('Has RT'): str,
     Optional('RT units'): str
     }
    )

ms_error_message = partial(create_error_message, 'MS')
ms_schema = Schema(
    {
        "INSTRUMENT_NAME": str,
        "INSTRUMENT_TYPE": str,
        "MS_TYPE": str,
        "ION_MODE": str,
        Optional("CAPILLARY_TEMPERATURE"): And(str, UnitChecker(['°C', 'C']), True, ms_error_message('CAPILLARY_TEMPERATURE', True)),
        Optional("CAPILLARY_VOLTAGE"): And(str, UnitChecker(['V', 'kV'], error_message = ms_error_message('CAPILLARY_VOLTAGE'))),
        Optional("COLLISION_ENERGY"): str,
        Optional("COLLISION_GAS"): And(str, SetChecker(['Nitrogen', 'Argon'])),
        Optional("DRY_GAS_FLOW"): And(str, UnitChecker(['L/hr', 'L/min'], error_message = ms_error_message('DRY_GAS_FLOW'))),
        Optional("DRY_GAS_TEMP"): And(str, UnitChecker(['°C', 'C'], error_message = ms_error_message('DRY_GAS_TEMP'))),
        Optional("FRAGMENT_VOLTAGE"): And(str, UnitChecker(['V'], error_message = ms_error_message('FRAGMENT_VOLTAGE'))),
        Optional("FRAGMENTATION_METHOD"): str,
        Optional("GAS_PRESSURE"): And(str, UnitChecker(['psi', 'psig', 'bar', 'kPa'], error_message = ms_error_message('GAS_PRESSURE'))),
        Optional("HELIUM_FLOW"): And(str, UnitChecker(['mL/min'], error_message = ms_error_message('HELIUM_FLOW'))),
        Optional("ION_SOURCE_TEMPERATURE"): And(str, UnitChecker(['°C', 'C'], error_message = ms_error_message('ION_SOURCE_TEMPERATURE'))),
        Optional("ION_SPRAY_VOLTAGE"): And(str, UnitChecker(['V', 'kV'], error_message = ms_error_message('ION_SPRAY_VOLTAGE'))),
        Optional("IONIZATION"): str,
        Optional("IONIZATION_ENERGY"): And(str, UnitChecker(['eV'], error_message = ms_error_message('IONIZATION_ENERGY'))),
        Optional("IONIZATION_POTENTIAL"): str,
        Optional("MASS_ACCURACY"): str,
        Optional("PRECURSOR_TYPE"): str,
        Optional("REAGENT_GAS"): str,
        Optional("SOURCE_TEMPERATURE"): And(str, UnitChecker(['°C', 'C'], error_message = ms_error_message('SOURCE_TEMPERATURE'))),
        Optional("SPRAY_VOLTAGE"): And(str, UnitChecker(['kV'], error_message = ms_error_message('SPRAY_VOLTAGE'))),
        Optional("ACTIVATION_PARAMETER"): str,
        Optional("ACTIVATION_TIME"): And(str, UnitChecker(['ms'], error_message = ms_error_message('ACTIVATION_TIME'))),
        Optional("ATOM_GUN_CURRENT"): str,
        Optional("AUTOMATIC_GAIN_CONTROL"): str,
        Optional("BOMBARDMENT"): str,
        Optional("CDL_SIDE_OCTOPOLES_BIAS_VOLTAGE"): str,
        Optional("CDL_TEMPERATURE"): str,
        Optional("DATAFORMAT"): str,
        Optional("DESOLVATION_GAS_FLOW"): And(str, UnitChecker(['L/hr', 'L/min'], error_message = ms_error_message('DESOLVATION_GAS_FLOW'))),
        Optional("DESOLVATION_TEMPERATURE"): And(str, UnitChecker(['°C', 'C'], error_message = ms_error_message('DESOLVATION_TEMPERATURE'))),
        Optional("INTERFACE_VOLTAGE"): str,
        Optional("IT_SIDE_OCTOPOLES_BIAS_VOLTAGE"): str,
        Optional("LASER"): str,
        Optional("MATRIX"): str,
        Optional("NEBULIZER"): str,
        Optional("OCTPOLE_VOLTAGE"): And(str, UnitChecker(['V'], error_message = ms_error_message('OCTPOLE_VOLTAGE'))),
        Optional("PROBE_TIP"): str,
        Optional("RESOLUTION_SETTING"): str,
        Optional("SAMPLE_DRIPPING"): str,
        Optional("SCAN_RANGE_MOVERZ"): str,
        Optional("SCANNING"): str,
        Optional("SCANNING_CYCLE"): str,
        Optional("SCANNING_RANGE"): str,
        Optional("SKIMMER_VOLTAGE"): And(str, UnitChecker(['V'], error_message = ms_error_message('SKIMMER_VOLTAGE'))),
        Optional("TUBE_LENS_VOLTAGE"): str,
        Optional("MS_COMMENTS"): str,  # changed to optional mwTab File Format Spec. 1.5
        Optional('MS_RESULTS_FILE'): results_file_schema
    }
)

nmr_error_message = partial(create_error_message, 'NM')
nmr_schema = Schema(
    {
        "INSTRUMENT_NAME": str,
        "INSTRUMENT_TYPE": str,
        "NMR_EXPERIMENT_TYPE": str,
        Optional("NMR_COMMENTS"): str,
        Optional("FIELD_FREQUENCY_LOCK"): str,
        Optional("STANDARD_CONCENTRATION"): And(str, UnitChecker(['mM'], error_message = nmr_error_message('STANDARD_CONCENTRATION'))),
        "SPECTROMETER_FREQUENCY": And(str, UnitChecker(['MHz'], error_message = nmr_error_message('SPECTROMETER_FREQUENCY'))),
        Optional("NMR_PROBE"): str,
        Optional("NMR_SOLVENT"): str,
        Optional("NMR_TUBE_SIZE"): str,
        Optional("SHIMMING_METHOD"): str,
        Optional("PULSE_SEQUENCE"): str,
        Optional("WATER_SUPPRESSION"): str,
        Optional("PULSE_WIDTH"): str,
        # TODO ask Hunter if dB is okay unit here.
        Optional("POWER_LEVEL"): And(str, UnitChecker(['W'], error_message = nmr_error_message('POWER_LEVEL'))),
        Optional("RECEIVER_GAIN"): And(str, NumChecker(nmr_error_message('RECEIVER_GAIN', False, True))),
        Optional("OFFSET_FREQUENCY"): And(str, UnitChecker(['ppm', 'Hz'], error_message = nmr_error_message('OFFSET_FREQUENCY'))),
        # TODO ask about dB.
        Optional("PRESATURATION_POWER_LEVEL"): And(str, UnitChecker(['W'], error_message = nmr_error_message('PRESATURATION_POWER_LEVEL'))),
        Optional("CHEMICAL_SHIFT_REF_CPD"): str,
        Optional("TEMPERATURE"): And(str, UnitChecker(['°C', 'C', 'K'], error_message = nmr_error_message('TEMPERATURE'))),
        Optional("NUMBER_OF_SCANS"): And(str, NumChecker(nmr_error_message('NUMBER_OF_SCANS', False, True))),
        Optional("DUMMY_SCANS"): str,
        Optional("ACQUISITION_TIME"): And(str, UnitChecker(['s'], error_message = nmr_error_message('ACQUISITION_TIME'))),
        Optional("RELAXATION_DELAY"): And(str, UnitChecker(['s', 'ms', 'us', 'μs'], error_message = nmr_error_message('RELAXATION_DELAY'))),
        Optional("SPECTRAL_WIDTH"): And(str, UnitChecker(['ppm', 'Hz'], error_message = nmr_error_message('SPECTRAL_WIDTH'))),
        Optional("NUM_DATA_POINTS_ACQUIRED"): And(str, NumChecker(nmr_error_message('NUM_DATA_POINTS_ACQUIRED', False, True))),
        Optional("REAL_DATA_POINTS"): str,
        Optional("LINE_BROADENING"): And(str, UnitChecker(['Hz'], error_message = nmr_error_message('LINE_BROADENING'))),
        Optional("ZERO_FILLING"): str,
        Optional("APODIZATION"): str,
        Optional("BASELINE_CORRECTION_METHOD"): str,
        Optional("CHEMICAL_SHIFT_REF_STD"): str,
        # TODO ask Hunter if NA is okay here. This seems like it could be more acceptable here.
        Optional("BINNED_INCREMENT"): And(str, UnitChecker(['ppm'], error_message = nmr_error_message('BINNED_INCREMENT'))),
        Optional("BINNED_DATA_NORMALIZATION_METHOD"): str,
        Optional("BINNED_DATA_PROTOCOL_FILE"): str,
        Optional("BINNED_DATA_CHEMICAL_SHIFT_RANGE"): str,
        Optional("BINNED_DATA_EXCLUDED_RANGE"): str,
        Optional("NMR_RESULTS_FILE"): results_file_schema
    }
)

data_schema = Schema(
    [
        {
            Or("Metabolite", "Bin range(ppm)", only_one=True): str,
            Optional(str): str,
        },
    ]
)

extended_schema = Schema(
    [
        {
            "Metabolite": str,
            Optional(str): str,
            "sample_id": str
        },
    ]
)

ms_metabolite_data_schema = Schema(
    {
        "Units": str,
        "Data": data_schema,
        Optional("Metabolites"): data_schema,
        Optional("Extended"): extended_schema
    }
)

nmr_binned_data_schema = Schema(
    {
        "Units": str,
        "Data": data_schema
    }
)

section_schema_mapping = {
    "METABOLOMICS WORKBENCH": metabolomics_workbench_schema,
    "PROJECT": project_schema,
    "STUDY": study_schema,
    "SUBJECT": subject_schema,
    "SUBJECT_SAMPLE_FACTORS": subject_sample_factors_schema,
    "COLLECTION": collection_schema,
    "TREATMENT": treatment_schema,
    "SAMPLEPREP": sampleprep_schema,
    "CHROMATOGRAPHY": chromatography_schema,
    "ANALYSIS": analysis_schema,
    "MS": ms_schema,
    "NM": nmr_schema,
    "MS_METABOLITE_DATA": ms_metabolite_data_schema,
    "NMR_METABOLITE_DATA": ms_metabolite_data_schema,
    "NMR_BINNED_DATA": nmr_binned_data_schema,
}


base_schema = Schema({
 "METABOLOMICS WORKBENCH": Or(dict, OrderedDict),
 "PROJECT": Or(dict, OrderedDict),
 "STUDY": Or(dict, OrderedDict),
 "SUBJECT": Or(dict, OrderedDict),
 "SUBJECT_SAMPLE_FACTORS": list,
 "COLLECTION": Or(dict, OrderedDict),
 "TREATMENT": Or(dict, OrderedDict),
 "SAMPLEPREP": Or(dict, OrderedDict),
 "ANALYSIS": Or(dict, OrderedDict),
 Optional("CHROMATOGRAPHY"): Or(dict, OrderedDict),
 Optional("MS_METABOLITE_DATA"): Or(dict, OrderedDict),
 Optional("NMR_METABOLITE_DATA"): Or(dict, OrderedDict),
 Optional("NMR_BINNED_DATA"): Or(dict, OrderedDict),
 Or("MS", "NM", only_one=True): Or(dict, OrderedDict),
 })

# _ms_base = deepcopy(_common_base)
# _ms_base["MS"] = Or(dict, OrderedDict)
# _ms_base["MS_METABOLITE_DATA"] = Or(dict, OrderedDict)
# _ms_base["CHROMATOGRAPHY"] = Or(dict, OrderedDict)

# _nmr_base = deepcopy(_common_base)
# _nmr_base["NM"] = Or(dict, OrderedDict)
# _nmr_base["NMR_METABOLITE_DATA"] = Or(dict, OrderedDict)

# _nmr_binned_base = deepcopy(_common_base)
# _nmr_binned_base["NM"] = Or(dict, OrderedDict)
# _nmr_binned_base["NMR_BINNED_DATA"] = Or(dict, OrderedDict)

# base_schema = Or(
#     Schema(_ms_base),
#     Schema(_nmr_base),
#     Schema(_nmr_binned_base),
#     )


