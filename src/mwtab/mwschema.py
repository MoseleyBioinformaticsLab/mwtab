#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
mwtab.mwschema
~~~~~~~~~~~~~~

This module provides schema definitions for different sections of the
``mwTab`` Metabolomics Workbench format.
"""

from copy import deepcopy
from functools import partial

import jsonschema

from . import metadata_column_matching



NA_VALUES = ['', '-', '−', '--', '---', 
             'NA', 'na', 'n.a.', 'N.A.', 'n/a', 'N/A', '#N/A', 'NaN', 'nan', 
             'null', 'Null', 'NULL', 'none', 'None',
             'unspecified', 'Unspecified']


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

def create_unit_error_message(validation_error: jsonschema.exceptions.ValidationError, 
                              mwtab_format: str = 'json',
                              can_be_range: bool = False, 
                              no_units: bool = False, 
                              units: list[str]|None = None) -> str:
    """Generate the error message for mwTab subsections that fail the unit regex.
    
    The idea for this function is that you include it in the jsonschema under a keyword that 
    is not reserved by jsonschema, and then your custom error handler will look for that 
    keyword and execute whatever function it finds there, passing the error object in as 
    the first parameter. You can combine this with the "partial" function from functools 
    to set the other parameters of this function as needed.
    
    Args:
        validation_error: The ValidationError created by jsonschema. Used to fill in some text in the message.
        mwtab_format: A string that must be either "mwtab" or "json". Customizes the error message for the format.
        can_be_range: If True, the regular expression used to validate could have matched a number range, so the message is modified to note that.
        no_units: If True, the regular expression did not require units to be present, so the message is modified to note that.
        units: If the regular expression required units, pass them in with this parameter so the message will indicate the allowed units.
    
    Returns:
        A completed string error message.
    """
    section = validation_error.relative_path[-2]
    subsection = validation_error.relative_path[-1]
    
    if mwtab_format == 'mwtab':
        format_string = f'for the subsection, "{subsection}", in the "{section}" section'
    else:
        format_string = f'in ["{section}"]["{subsection}"]'
    
    if can_be_range:
        range_string = ' or range (ex. "5-6") '
    else:
        range_string = ''
    
    if no_units:
        unit_string = '.'
    else:
        unit_string = f'followed by a space with a unit (ex. "5 V") from the following list: {units}.'
    
    message = (f'The value, "{validation_error.instance}", {format_string} '
               f'should be a {"unitless " if no_units else ""}number{range_string}{unit_string}')
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

def _create_num_regex_and_message_func(can_be_range: bool = False) -> str:
    """Simple wrapper for DRY purposes.
    
    Creating the regular expression and validation error message in 1 function mixes too many 
    concerns into 1 function, so they are split into 2 and this function serves as a convenience 
    to pass them into a jsonschema easily. To this end the return is in dictionary form with the 
    intention for it to be unpacked into the jsonschema. For example, 
    {**_create_num_regex_and_message_func(True)}
    
    Args:
        can_be_range: If true, the regex will match a number range and the message will be slightly different.
    
    Returns:
        A dicitonary {'pattern': regex, 'message_func': message_function}.
    """
    regex = '^\d+$'
    message_function = partial(create_unit_error_message, can_be_range = can_be_range, no_units = True, units = None)
    return {'pattern': regex, 'message_func': message_function}


def create_ID_error_message(validation_error: jsonschema.exceptions.ValidationError, 
                            mwtab_format: str = 'json') -> str:
    """Generate the error message for IDs in the mwtab header.
    
    Similar to create_unit_error_message.
    
    Args:
        validation_error: The ValidationError created by jsonschema. Used to fill in some text in the message.
        mwtab_format: A string that must be either "mwtab" or "json". Customizes the error message for the format.
    
    Returns:
        A completed string error message.
    """
    section = validation_error.relative_path[-2]
    subsection = validation_error.relative_path[-1]
    
    if subsection == 'STUDY_ID':
        abbrev = 'ST'
    elif subsection == 'PROJECT_ID':
        abbrev = 'PR'
    elif subsection == 'ANALYSIS_ID':
        abbrev = 'AN'
    else:
        raise ValueError('Unknown subsection.')
    
    if mwtab_format == 'mwtab':
        format_string = f'for "{subsection}" in the file header'
    else:
        format_string = f'in ["{section}"]["{subsection}"]'
    
    message = (f'The value, "{validation_error.instance}", {format_string} '
               f'must be the letters "{abbrev}" followed by 6 numbers. Ex. "{abbrev}001405".')
    return message

def create_numeric_error_message(validation_error: jsonschema.exceptions.ValidationError, 
                                 mwtab_format: str = 'json') -> str:
    """Generate the error message for numbers in the mwtab header.
    
    Similar to create_unit_error_message.
    
    Args:
        validation_error: The ValidationError created by jsonschema. Used to fill in some text in the message.
        mwtab_format: A string that must be either "mwtab" or "json". Customizes the error message for the format.
    
    Returns:
        A completed string error message.
    """
    section = validation_error.relative_path[-2]
    subsection = validation_error.relative_path[-1]
    
    if mwtab_format == 'mwtab':
        format_string = f'for "{subsection}" in the file header'
    else:
        format_string = f'in ["{section}"]["{subsection}"]'
    
    message = (f'The value, "{validation_error.instance}", {format_string} '
               'must be a positive integer. Ex. "1" or "2051".')
    return message







metabolomics_workbench_schema = \
{'type': 'object',
 'properties': {'STUDY_ID': {'type': 'string', 'pattern': r'^ST\d{6}$', 'message_func': create_ID_error_message},
                'ANALYSIS_ID': {'type': 'string', 'pattern': r'^AN\d{6}$', 'message_func': create_ID_error_message},
                'VERSION': {'type': 'string', 'pattern': r'^\d+$', 'message_func': create_numeric_error_message},
                'CREATED_ON': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'PROJECT_ID': {'type': 'string', 'pattern': r'^PR\d{6}$', 'message_func': create_ID_error_message},
                'HEADER': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'DATATRACK_ID': {'type': 'string', 'pattern': r'^\d+$', 'message_func': create_numeric_error_message},
                'filename': {'type': 'string'}},
 'required': ['VERSION', 'CREATED_ON'],
 'additionalProperties': False}

project_schema = \
{'type': 'object',
 'properties': {'PROJECT_TITLE': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'PROJECT_TYPE': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'PROJECT_SUMMARY': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'INSTITUTE': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'DEPARTMENT': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'LABORATORY': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'LAST_NAME': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'FIRST_NAME': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'ADDRESS': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'EMAIL': {'type': 'string', 'format': 'email'},
                'PHONE': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'FUNDING_SOURCE': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'PROJECT_COMMENTS': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'PUBLICATIONS': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'CONTRIBUTORS': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                # TODO test this across the files. I think some are going to have "doi.org" in front.
                'DOI': {'type': 'string', 'pattern': r'10\.\d{4,9}/[-._;()/:a-z0-9A-Z]+', 'error_message': 'The "DOI" in the "PROJECT" section does not appear to be valid.'}},
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
 'properties': {'STUDY_TITLE': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'STUDY_TYPE': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'STUDY_SUMMARY': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'INSTITUTE': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'DEPARTMENT': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'LABORATORY': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'LAST_NAME': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'FIRST_NAME': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'ADDRESS': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'EMAIL': {'type': 'string', 'format': 'email'},
                'PHONE': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'SUBMIT_DATE': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'NUM_GROUPS': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'TOTAL_SUBJECTS': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'NUM_MALES': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'NUM_FEMALES': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'STUDY_COMMENTS': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'PUBLICATIONS': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}}},
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
 'properties': {'SUBJECT_TYPE': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'SUBJECT_SPECIES': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'TAXONOMY_ID': {'type': 'string', 'pattern': metadata_column_matching.make_list_regex(r'\d+', r'(,|;|\||/)'), 'error_message': '"TAXONOMY_ID" in the "SUBJECT" section must be a number or list of numbers.'},
                'GENOTYPE_STRAIN': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'AGE_OR_AGE_RANGE': {'type': 'string', **_create_unit_regex_and_message_func(['weeks', 'days', 'months', 'years'], True)},
                'WEIGHT_OR_WEIGHT_RANGE': {'type': 'string', **_create_unit_regex_and_message_func(['g', 'mg', 'kg', 'lbs'], True)},
                'HEIGHT_OR_HEIGHT_RANGE': {'type': 'string', **_create_unit_regex_and_message_func(['cm', 'in'], True)},
                'GENDER': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'HUMAN_RACE': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'HUMAN_ETHNICITY': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'HUMAN_TRIAL_TYPE': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'HUMAN_LIFESTYLE_FACTORS': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'HUMAN_MEDICATIONS': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'HUMAN_PRESCRIPTION_OTC': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'HUMAN_SMOKING_STATUS': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'HUMAN_ALCOHOL_DRUG_USE': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'HUMAN_NUTRITION': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'HUMAN_INCLUSION_CRITERIA': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'HUMAN_EXCLUSION_CRITERIA': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'ANIMAL_ANIMAL_SUPPLIER': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'ANIMAL_HOUSING': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'ANIMAL_LIGHT_CYCLE': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'ANIMAL_FEED': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'ANIMAL_WATER': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'ANIMAL_INCLUSION_CRITERIA': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'CELL_BIOSOURCE_OR_SUPPLIER': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'CELL_STRAIN_DETAILS': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'SUBJECT_COMMENTS': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'CELL_PRIMARY_IMMORTALIZED': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'CELL_PASSAGE_NUMBER': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'CELL_COUNTS': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'SPECIES_GROUP': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}}},
 'required': ['SUBJECT_TYPE', 'SUBJECT_SPECIES'],
 'additionalProperties': False}

subject_sample_factors_schema = \
{'type': 'array',
 'items': {'type': 'object',
           'properties': {'Subject ID': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                          'Sample ID': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                          'Factors': {'type': 'object'},
                          'Additional sample data': {'type': 'object',
                                                     'properties': {'RAW_FILE_NAME': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}}},
                                                     'additionalProperties': True}},
           'required': ['Subject ID', 'Sample ID', 'Factors'],
           'additionalProperties': False}}    

collection_schema = \
{'type': 'object',
 'properties': {'COLLECTION_SUMMARY': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'COLLECTION_PROTOCOL_ID': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'COLLECTION_PROTOCOL_FILENAME': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'COLLECTION_PROTOCOL_COMMENTS': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'SAMPLE_TYPE': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'COLLECTION_METHOD': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'COLLECTION_LOCATION': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'COLLECTION_FREQUENCY': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'COLLECTION_DURATION': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'COLLECTION_TIME': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'VOLUMEORAMOUNT_COLLECTED': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'STORAGE_CONDITIONS': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'COLLECTION_VIALS': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'STORAGE_VIALS': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'COLLECTION_TUBE_TEMP': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'ADDITIVES': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                # TODO Ask Hunter if these should be lowercase or leave them as is.
                'BLOOD_SERUM_OR_PLASMA': {'type': 'string', 'enum':['Blood', 'Serum', 'Plasma']},
                'TISSUE_CELL_IDENTIFICATION': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'TISSUE_CELL_QUANTITY_TAKEN': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}}},
 'required': ['COLLECTION_SUMMARY'],
 'additionalProperties': False}

treatment_schema = \
{'type': 'object',
 'properties': {'TREATMENT_SUMMARY': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'TREATMENT_PROTOCOL_ID': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'TREATMENT_PROTOCOL_FILENAME': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'TREATMENT_PROTOCOL_COMMENTS': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'TREATMENT': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'TREATMENT_COMPOUND': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'TREATMENT_ROUTE': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'TREATMENT_DOSE': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'TREATMENT_DOSEVOLUME': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'TREATMENT_DOSEDURATION': {'type': 'string', **_create_unit_regex_and_message_func(['h', 'weeks', 'days'], True)},
                'TREATMENT_VEHICLE': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'ANIMAL_VET_TREATMENTS': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'ANIMAL_ANESTHESIA': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'ANIMAL_ACCLIMATION_DURATION': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'ANIMAL_FASTING': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'ANIMAL_ENDP_EUTHANASIA': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'ANIMAL_ENDP_TISSUE_COLL_LIST': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'ANIMAL_ENDP_TISSUE_PROC_METHOD': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'ANIMAL_ENDP_CLINICAL_SIGNS': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'HUMAN_FASTING': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'HUMAN_ENDP_CLINICAL_SIGNS': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'CELL_STORAGE': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'CELL_GROWTH_CONTAINER': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'CELL_GROWTH_CONFIG': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'CELL_GROWTH_RATE': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'CELL_INOC_PROC': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'CELL_MEDIA': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'CELL_ENVIR_COND': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'CELL_HARVESTING': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'PLANT_GROWTH_SUPPORT': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'PLANT_GROWTH_LOCATION': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'PLANT_PLOT_DESIGN': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'PLANT_LIGHT_PERIOD': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'PLANT_HUMIDITY': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'PLANT_TEMP': {'type': 'string', **_create_unit_regex_and_message_func(['°C', 'C'], True)},
                'PLANT_WATERING_REGIME': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'PLANT_NUTRITIONAL_REGIME': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'PLANT_ESTAB_DATE': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'PLANT_HARVEST_DATE': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'PLANT_GROWTH_STAGE': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'PLANT_METAB_QUENCH_METHOD': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'PLANT_HARVEST_METHOD': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'PLANT_STORAGE': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'CELL_PCT_CONFLUENCE': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'CELL_MEDIA_LASTCHANGED': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}}},
 'required': ['TREATMENT_SUMMARY'],
 'additionalProperties': False}

sampleprep_schema = \
{'type': 'object',
 'properties': {'SAMPLEPREP_SUMMARY': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'SAMPLEPREP_PROTOCOL_ID': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'SAMPLEPREP_PROTOCOL_FILENAME': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'SAMPLEPREP_PROTOCOL_COMMENTS': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'PROCESSING_METHOD': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'PROCESSING_STORAGE_CONDITIONS': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'EXTRACTION_METHOD': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'EXTRACT_CONCENTRATION_DILUTION': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'EXTRACT_ENRICHMENT': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'EXTRACT_CLEANUP': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'EXTRACT_STORAGE': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'SAMPLE_RESUSPENSION': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'SAMPLE_DERIVATIZATION': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'SAMPLE_SPIKING': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'ORGAN': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'ORGAN_SPECIFICATION': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'CELL_TYPE': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'SUBCELLULAR_LOCATION': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}}},
 'required': ['SAMPLEPREP_SUMMARY'],
 'additionalProperties': False}

chromatography_schema = \
{'type': 'object',
 'properties': {'CHROMATOGRAPHY_SUMMARY': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'CHROMATOGRAPHY_TYPE': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'INSTRUMENT_NAME': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'COLUMN_NAME': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'FLOW_GRADIENT': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'FLOW_RATE': {'type': 'string', **_create_unit_regex_and_message_func(['mL/min', 'uL/min', 'μL/min'], True)},
                'COLUMN_TEMPERATURE': {'type': 'string', **_create_unit_regex_and_message_func(['°C', 'C'], True)},
                'METHODS_FILENAME': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'SAMPLE_INJECTION': {'type': 'string', **_create_unit_regex_and_message_func(['μL', 'uL'])},
                'SOLVENT_A': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'SOLVENT_B': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'METHODS_ID': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'COLUMN_PRESSURE': {'type': 'string', **_create_unit_regex_and_message_func(['psi', 'bar'], True)},
                # TODO ask Hunter about a special case for "room temperature".
                'INJECTION_TEMPERATURE': {'type': 'string', **_create_unit_regex_and_message_func(['°C', 'C'], True)},
                'INTERNAL_STANDARD': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'INTERNAL_STANDARD_MT': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'RETENTION_INDEX': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'RETENTION_TIME': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'SAMPLING_CONE': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'ANALYTICAL_TIME': {'type': 'string', **_create_unit_regex_and_message_func(['min'], True)},
                'CAPILLARY_VOLTAGE': {'type': 'string', **_create_unit_regex_and_message_func(['V', 'kV'])},
                'MIGRATION_TIME': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'OVEN_TEMPERATURE': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'PRECONDITIONING': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'RUNNING_BUFFER': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'RUNNING_VOLTAGE': {'type': 'string', **_create_unit_regex_and_message_func(['V', 'kV'])},
                'SHEATH_LIQUID': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'TIME_PROGRAM': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'TRANSFERLINE_TEMPERATURE': {'type': 'string', **_create_unit_regex_and_message_func(['°C', 'C'])},
                'WASHING_BUFFER': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'WEAK_WASH_SOLVENT_NAME': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'WEAK_WASH_VOLUME': {'type': 'string', **_create_unit_regex_and_message_func(['μL', 'uL'])},
                'STRONG_WASH_SOLVENT_NAME': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'STRONG_WASH_VOLUME': {'type': 'string', **_create_unit_regex_and_message_func(['μL', 'uL'])},
                'TARGET_SAMPLE_TEMPERATURE': {'type': 'string', **_create_unit_regex_and_message_func(['°C', 'C'])},
                'SAMPLE_LOOP_SIZE': {'type': 'string', **_create_unit_regex_and_message_func(['μL', 'uL'])},
                'SAMPLE_SYRINGE_SIZE': {'type': 'string', **_create_unit_regex_and_message_func(['μL', 'uL'])},
                'RANDOMIZATION_ORDER': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'CHROMATOGRAPHY_COMMENTS': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}}},
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
                'LABORATORY_NAME': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'ACQUISITION_DATE': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'SOFTWARE_VERSION': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'OPERATOR_NAME': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'DETECTOR_TYPE': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'ANALYSIS_PROTOCOL_FILE': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'ACQUISITION_PARAMETERS_FILE': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'PROCESSING_PARAMETERS_FILE': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'DATA_FORMAT': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'ACQUISITION_ID': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'ACQUISITION_TIME': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'ANALYSIS_COMMENTS': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'ANALYSIS_DISPLAY': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'INSTRUMENT_NAME': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'INSTRUMENT_PARAMETERS_FILE': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'NUM_FACTORS': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'NUM_METABOLITES': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'PROCESSED_FILE': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'RANDOMIZATION_ORDER': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'RAW_FILE': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}}},
 'required': ['ANALYSIS_TYPE'],
 'additionalProperties': False}

results_file_schema = \
{'type': 'object',
 'properties': {'filename': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
  'UNITS': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
  'Has m/z': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
  'Has RT': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
  'RT units': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}}},
 'required': ['filename', 'UNITS'],
 'additionalProperties': False}

ms_schema = \
{'type': 'object',
 'properties': {'INSTRUMENT_NAME': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'INSTRUMENT_TYPE': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'MS_TYPE': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'ION_MODE': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'CAPILLARY_TEMPERATURE': {'type': 'string', **_create_unit_regex_and_message_func(['°C', 'C'], True)},
                'CAPILLARY_VOLTAGE': {'type': 'string', **_create_unit_regex_and_message_func(['V', 'kV'])},
                'COLLISION_ENERGY': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'COLLISION_GAS': {'type': 'string', 'enum': ['Nitrogen', 'Argon']},
                'DRY_GAS_FLOW': {'type': 'string', **_create_unit_regex_and_message_func(['L/hr', 'L/min'])},
                'DRY_GAS_TEMP': {'type': 'string', **_create_unit_regex_and_message_func(['°C', 'C'])},
                'FRAGMENT_VOLTAGE': {'type': 'string', **_create_unit_regex_and_message_func(['V'])},
                'FRAGMENTATION_METHOD': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'GAS_PRESSURE': {'type': 'string', **_create_unit_regex_and_message_func(['psi', 'psig', 'bar', 'kPa'])},
                'HELIUM_FLOW': {'type': 'string', **_create_unit_regex_and_message_func(['mL/min'])},
                'ION_SOURCE_TEMPERATURE': {'type': 'string', **_create_unit_regex_and_message_func(['°C', 'C'])},
                'ION_SPRAY_VOLTAGE': {'type': 'string', **_create_unit_regex_and_message_func(['V', 'kV'])},
                'IONIZATION': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'IONIZATION_ENERGY': {'type': 'string', **_create_unit_regex_and_message_func(['eV'])},
                'IONIZATION_POTENTIAL': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'MASS_ACCURACY': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'PRECURSOR_TYPE': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'REAGENT_GAS': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'SOURCE_TEMPERATURE': {'type': 'string', **_create_unit_regex_and_message_func(['°C', 'C'])},
                'SPRAY_VOLTAGE': {'type': 'string', **_create_unit_regex_and_message_func(['kV'])},
                'ACTIVATION_PARAMETER': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'ACTIVATION_TIME': {'type': 'string', **_create_unit_regex_and_message_func(['ms'])},
                'ATOM_GUN_CURRENT': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'AUTOMATIC_GAIN_CONTROL': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'BOMBARDMENT': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'CDL_SIDE_OCTOPOLES_BIAS_VOLTAGE': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'CDL_TEMPERATURE': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'DATAFORMAT': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'DESOLVATION_GAS_FLOW': {'type': 'string', **_create_unit_regex_and_message_func(['L/hr', 'L/min'])},
                'DESOLVATION_TEMPERATURE': {'type': 'string', **_create_unit_regex_and_message_func(['°C', 'C'])},
                'INTERFACE_VOLTAGE': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'IT_SIDE_OCTOPOLES_BIAS_VOLTAGE': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'LASER': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'MATRIX': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'NEBULIZER': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'OCTPOLE_VOLTAGE': {'type': 'string', **_create_unit_regex_and_message_func(['V'])},
                'PROBE_TIP': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'RESOLUTION_SETTING': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'SAMPLE_DRIPPING': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'SCAN_RANGE_MOVERZ': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'SCANNING': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'SCANNING_CYCLE': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'SCANNING_RANGE': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'SKIMMER_VOLTAGE': {'type': 'string', **_create_unit_regex_and_message_func(['V'])},
                'TUBE_LENS_VOLTAGE': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'MS_COMMENTS': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'MS_RESULTS_FILE': results_file_schema},
 'required': ['INSTRUMENT_NAME', 'INSTRUMENT_TYPE', 'MS_TYPE', 'ION_MODE'],
 'additionalProperties': False}

nmr_schema = \
{'type': 'object',
 'properties': {'INSTRUMENT_NAME': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'INSTRUMENT_TYPE': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'NMR_EXPERIMENT_TYPE': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'NMR_COMMENTS': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'FIELD_FREQUENCY_LOCK': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'STANDARD_CONCENTRATION': {'type': 'string', **_create_unit_regex_and_message_func(['mM'])},
                'SPECTROMETER_FREQUENCY': {'type': 'string', **_create_unit_regex_and_message_func(['MHz'])},
                'NMR_PROBE': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'NMR_SOLVENT': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'NMR_TUBE_SIZE': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'SHIMMING_METHOD': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'PULSE_SEQUENCE': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'WATER_SUPPRESSION': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'PULSE_WIDTH': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                # TODO ask Hunter if dB is okay unit here.
                'POWER_LEVEL': {'type': 'string', **_create_unit_regex_and_message_func(['W'])},
                'RECEIVER_GAIN': {'type': 'string', **_create_num_regex_and_message_func(False)},
                'OFFSET_FREQUENCY': {'type': 'string', **_create_unit_regex_and_message_func(['ppm', 'Hz'])},
                # TODO ask about dB.
                'PRESATURATION_POWER_LEVEL': {'type': 'string', **_create_unit_regex_and_message_func(['W'])},
                'CHEMICAL_SHIFT_REF_CPD': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'TEMPERATURE': {'type': 'string', **_create_unit_regex_and_message_func(['°C', 'C', 'K'])},
                'NUMBER_OF_SCANS': {'type': 'string', **_create_num_regex_and_message_func(False)},
                'DUMMY_SCANS': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'ACQUISITION_TIME': {'type': 'string', **_create_unit_regex_and_message_func(['s'])},
                'RELAXATION_DELAY': {'type': 'string', **_create_unit_regex_and_message_func(['s', 'ms', 'us', 'μs'])},
                'SPECTRAL_WIDTH': {'type': 'string', **_create_unit_regex_and_message_func(['ppm', 'Hz'])},
                'NUM_DATA_POINTS_ACQUIRED': {'type': 'string', **_create_num_regex_and_message_func(False)},
                'REAL_DATA_POINTS': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'LINE_BROADENING': {'type': 'string', **_create_unit_regex_and_message_func(['Hz'])},
                'ZERO_FILLING': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'APODIZATION': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'BASELINE_CORRECTION_METHOD': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'CHEMICAL_SHIFT_REF_STD': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'BINNED_INCREMENT': {'type': 'string', **_create_unit_regex_and_message_func(['ppm'])},
                'BINNED_DATA_NORMALIZATION_METHOD': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'BINNED_DATA_PROTOCOL_FILE': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'BINNED_DATA_CHEMICAL_SHIFT_RANGE': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'BINNED_DATA_EXCLUDED_RANGE': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'NMR_RESULTS_FILE': results_file_schema},
 'required': ['INSTRUMENT_NAME',
              'INSTRUMENT_TYPE',
              'NMR_EXPERIMENT_TYPE',
              'SPECTROMETER_FREQUENCY'],
 'additionalProperties': False}

# TODO consider removing the data sections because we validate these as tables, so checking each item in the array is 
# probably not necessary and might even lead to duplicate errors or errors that shouldn't be.
data_schema = \
{'type': 'array',
 'items': {'type': 'object',
           'properties': {'Metabolite': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                          'Bin range(ppm)': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}}},
           'required': [],
           'additionalProperties': True}}

extended_schema = \
{'type': 'array',
 'items': {'type': 'object',
           'properties': {'Metabolite': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                          'sample_id': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}}},
           'required': ['Metabolite', 'sample_id'],
           'additionalProperties': True}}

ms_metabolite_data_schema = \
{'type': 'object',
 'properties': {'Units': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
                'Data': data_schema,
                'Metabolites': data_schema,
                'Extended': extended_schema},
 'required': ['Units', 'Data'],
 'additionalProperties': False}

nmr_binned_data_schema = \
{'type': 'object',
 'properties': {'Units': {'type': 'string', 'minLength': 1, 'not':{'enum': NA_VALUES}},
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
    
    
