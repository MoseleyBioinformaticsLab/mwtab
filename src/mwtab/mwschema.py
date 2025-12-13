#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
mwtab.mwschema
~~~~~~~~~~~~~~

This module provides schema definitions for different sections of the
``mwTab`` Metabolomics Workbench format.
"""

from copy import deepcopy

from . import metadata_column_matching



NA_VALUES = ['', '-', '−', '--', '---', 
             'NA', 'na', 'n.a.', 'N.A.', 'n/a', 'N/A', '#N/A', 'NaN', 'nan', 
             'null', 'Null', 'NULL', 'none', 'None',
             'unspecified', 'Unspecified']
# 'NA' is a legitimate metabolite name that very rarely shows up.
METABOLITE_NA_VALUES = [value for value in NA_VALUES if value != 'NA']


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

def create_unit_error_message(can_be_range: bool = False, 
                              no_units: bool = False, 
                              units: list[str]|None = None,
                              integer: bool = False) -> str:
    """Generate the error message for mwTab subsections that fail the unit regex.
    
    Args:
        can_be_range: If True, the regular expression used to validate could have matched a number range, so the message is modified to note that.
        no_units: If True, the regular expression did not require units to be present, so the message is modified to note that.
        units: If the regular expression required units, pass them in with this parameter so the message will indicate the allowed units.
        integer: If True, only integers are allowed so the message will refer to intergers instead of numbers.
    
    Returns:
        A completed string error message.
    """
    if can_be_range:
        range_string = ' or range (ex. "5-6") '
    else:
        range_string = '' if no_units else ' '
    
    if no_units:
        unit_string = '.'
    else:
        unit_string = f'followed by a space with a unit (ex. "5 V") from the following list: {units}.'
    
    if integer:
        number = 'integer'
    else:
        number = 'number'
    
    message = (f' should be a {"unitless " if no_units else ""}{number}{range_string}{unit_string} '
               'Ignore this when more complicated descriptions are required.')
    return message

def _create_unit_regex_and_message(units: list[str], can_be_range: bool = False) -> str:
    """Simple wrapper for DRY purposes.
    
    Creating the regular expression and validation error message in 1 function mixes too many 
    concerns into 1 function, so they are split into 2 and this function serves as a convenience 
    to pass them into a jsonschema easily. To this end the return is in dictionary form with the 
    intention for it to be unpacked into the jsonschema. For example, 
    {**_create_unit_regex_and_message(['V'], True)}
    
    Args:
        units: A list of strings used to create the regex and error message.
        can_be_range: If true, the regex will match a number range and the message will be slightly different.
    
    Returns:
        A dicitonary {'pattern': regex, 'message_func': message_function}.
    """
    regex = create_units_regex(units, can_be_range)
    message = create_unit_error_message(can_be_range = can_be_range, no_units = False, units = units)
    return {'pattern': regex, 'pattern_custom_message': message}

def _create_num_regex_and_message(can_be_range: bool = False) -> str:
    """Simple wrapper for DRY purposes.
    
    Creating the regular expression and validation error message in 1 function mixes too many 
    concerns into 1 function, so they are split into 2 and this function serves as a convenience 
    to pass them into a jsonschema easily. To this end the return is in dictionary form with the 
    intention for it to be unpacked into the jsonschema. For example, 
    {**_create_num_regex_and_message(True)}
    
    Args:
        can_be_range: If true, the regex will match a number range and the message will be slightly different.
    
    Returns:
        A dicitonary {'pattern': regex, 'message_func': message_function}.
    """
    regex = r'^((\d+)|(\d*\.\d+))$'
    message = create_unit_error_message(can_be_range = can_be_range, no_units = True, units = None)
    return {'pattern': regex, 'pattern_custom_message': message}

def _create_int_regex_and_message(can_be_range: bool = False) -> str:
    """Simple wrapper for DRY purposes.
    
    Creating the regular expression and validation error message in 1 function mixes too many 
    concerns into 1 function, so they are split into 2 and this function serves as a convenience 
    to pass them into a jsonschema easily. To this end the return is in dictionary form with the 
    intention for it to be unpacked into the jsonschema. For example, 
    {**_create_num_regex_and_message(True)}
    
    Args:
        can_be_range: If true, the regex will match an integer range and the message will be slightly different.
    
    Returns:
        A dicitonary {'pattern': regex, 'message_func': message_function}.
    """
    regex = r'^\d+$'
    message = create_unit_error_message(can_be_range = can_be_range, no_units = True, units = None)
    return {'pattern': regex, 'pattern_custom_message': message}




ID_error_template = ' must be the letters "{abbrev}" followed by 6 numbers. Ex. "{abbrev}001405".'
numeric_error_template = ' must be a positive integer. Ex. "1" or "2051".'
metabolomics_workbench_schema = \
{'type': 'object',
 'properties': {'STUDY_ID': {'type': 'string', 'pattern': r'^ST\d{6}$', 'pattern_custom_message': ID_error_template.format(abbrev = 'ST')},
                'ANALYSIS_ID': {'type': 'string', 'pattern': r'^AN\d{6}$', 'pattern_custom_message': ID_error_template.format(abbrev = 'AN')},
                'VERSION': {'type': 'string', 'pattern': r'^\d+$', 'pattern_custom_message': numeric_error_template},
                'CREATED_ON': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'PROJECT_ID': {'type': 'string', 'pattern': r'^PR\d{6}$', 'pattern_custom_message': ID_error_template.format(abbrev = 'PR')},
                'HEADER': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'DATATRACK_ID': {'type': 'string', 'pattern': r'^\d+$', 'pattern_custom_message': numeric_error_template},
                'filename': {'type': 'string'}},
 'required': ['VERSION', 'CREATED_ON'],
 'additionalProperties': False}

project_schema = \
{'type': 'object',
 'properties': {'PROJECT_TITLE': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'PROJECT_TYPE': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'PROJECT_SUMMARY': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'INSTITUTE': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'DEPARTMENT': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'LABORATORY': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'LAST_NAME': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'FIRST_NAME': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'ADDRESS': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'EMAIL': {'type': 'string', 'format': 'email'},
                'PHONE': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'FUNDING_SOURCE': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'PROJECT_COMMENTS': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'PUBLICATIONS': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'CONTRIBUTORS': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'DOI': {'type': 'string', 'pattern': r'10\.\d{4,9}/[-._;()/:a-z0-9A-Z]+', 
                                          'pattern_custom_message': ' does not appear to be a valid DOI.'}},
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
 'properties': {'STUDY_TITLE': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'STUDY_TYPE': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'STUDY_SUMMARY': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'INSTITUTE': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'DEPARTMENT': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'LABORATORY': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'LAST_NAME': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'FIRST_NAME': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'ADDRESS': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'EMAIL': {'type': 'string', 'format': 'email'},
                'PHONE': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'SUBMIT_DATE': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'NUM_GROUPS': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'TOTAL_SUBJECTS': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'NUM_MALES': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'NUM_FEMALES': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'STUDY_COMMENTS': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'PUBLICATIONS': {'type': 'string', 'not':{'enum': NA_VALUES}}},
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
 'properties': {'SUBJECT_TYPE': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'SUBJECT_SPECIES': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'TAXONOMY_ID': {'type': 'string', 
                                'pattern': metadata_column_matching.make_list_regex(r'\d+', r'(,|;|\||/)', empty_string = True), 
                                'pattern_custom_message': ' must be a number or list of numbers.',
                                'not':{'enum': NA_VALUES}},
                'GENOTYPE_STRAIN': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'AGE_OR_AGE_RANGE': {'type': 'string', **_create_unit_regex_and_message(['weeks', 'days', 'months', 'years'], True)},
                'WEIGHT_OR_WEIGHT_RANGE': {'type': 'string', **_create_unit_regex_and_message(['g', 'mg', 'kg', 'lbs'], True)},
                'HEIGHT_OR_HEIGHT_RANGE': {'type': 'string', **_create_unit_regex_and_message(['cm', 'in'], True)},
                'GENDER': {'type': 'string', 'pattern': r'^((?i:male)|(?i:female)|(?i:male, female)|(?i:hermaphrodite)|N/A)$',
                           'pattern_custom_message': (' should be one of "Male", "Female", "Male, Female", "Hermaphrodite", or "N/A". '
                                                      'Ignore this when more complicated descriptions are required.')},
                'HUMAN_RACE': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'HUMAN_ETHNICITY': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'HUMAN_TRIAL_TYPE': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'HUMAN_LIFESTYLE_FACTORS': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'HUMAN_MEDICATIONS': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'HUMAN_PRESCRIPTION_OTC': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'HUMAN_SMOKING_STATUS': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'HUMAN_ALCOHOL_DRUG_USE': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'HUMAN_NUTRITION': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'HUMAN_INCLUSION_CRITERIA': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'HUMAN_EXCLUSION_CRITERIA': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'ANIMAL_ANIMAL_SUPPLIER': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'ANIMAL_HOUSING': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'ANIMAL_LIGHT_CYCLE': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'ANIMAL_FEED': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'ANIMAL_WATER': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'ANIMAL_INCLUSION_CRITERIA': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'CELL_BIOSOURCE_OR_SUPPLIER': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'CELL_STRAIN_DETAILS': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'SUBJECT_COMMENTS': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'CELL_PRIMARY_IMMORTALIZED': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'CELL_PASSAGE_NUMBER': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'CELL_COUNTS': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'SPECIES_GROUP': {'type': 'string', 'not':{'enum': NA_VALUES}}},
 'required': ['SUBJECT_TYPE', 'SUBJECT_SPECIES'],
 'additionalProperties': False}

subject_sample_factors_schema = \
{'type': 'array',
 'items': {'type': 'object',
           'properties': {'Subject ID': {'type': 'string'},
                          'Sample ID': {'type': 'string', 'not':{'enum': NA_VALUES}},
                          'Factors': {'type': 'object',
                                      'additionalProperties': {'type': 'string', 'minLength':1, 
                                                               'minLength_prefix':'Warning: ',
                                                               'minLength_custom_message': ' is missing a value.'}},
                          'Additional sample data': {'type': 'object',
                                                     'properties': {'RAW_FILE_NAME': {'type': 'string', 'not':{'enum': NA_VALUES}}},
                                                     'additionalProperties': {'type': 'string', 'minLength':1, 
                                                                              'minLength_prefix':'Warning: ',
                                                                              'minLength_custom_message': ' is missing a value.'}}},
           'required': ['Subject ID', 'Sample ID', 'Factors'],
           'additionalProperties': False}}    

collection_schema = \
{'type': 'object',
 'properties': {'COLLECTION_SUMMARY': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'COLLECTION_PROTOCOL_ID': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'COLLECTION_PROTOCOL_FILENAME': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'COLLECTION_PROTOCOL_COMMENTS': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'SAMPLE_TYPE': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'COLLECTION_METHOD': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'COLLECTION_LOCATION': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'COLLECTION_FREQUENCY': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'COLLECTION_DURATION': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'COLLECTION_TIME': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'VOLUMEORAMOUNT_COLLECTED': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'STORAGE_CONDITIONS': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'COLLECTION_VIALS': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'STORAGE_VIALS': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'COLLECTION_TUBE_TEMP': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'ADDITIVES': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'BLOOD_SERUM_OR_PLASMA': {'type': 'string', 'pattern': r'(?i)^(blood|plasma|serum)$',
                                          'pattern_custom_message': (' should be one of "Blood", "Plasma", or "Serum". '
                                                                     'Ignore this when more complicated descriptions are required.')},
                'TISSUE_CELL_IDENTIFICATION': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'TISSUE_CELL_QUANTITY_TAKEN': {'type': 'string', 'not':{'enum': NA_VALUES}}},
 'required': ['COLLECTION_SUMMARY'],
 'additionalProperties': False}

treatment_schema = \
{'type': 'object',
 'properties': {'TREATMENT_SUMMARY': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'TREATMENT_PROTOCOL_ID': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'TREATMENT_PROTOCOL_FILENAME': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'TREATMENT_PROTOCOL_COMMENTS': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'TREATMENT': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'TREATMENT_COMPOUND': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'TREATMENT_ROUTE': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'TREATMENT_DOSE': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'TREATMENT_DOSEVOLUME': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'TREATMENT_DOSEDURATION': {'type': 'string', **_create_unit_regex_and_message(['h', 'weeks', 'days'], True)},
                'TREATMENT_VEHICLE': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'ANIMAL_VET_TREATMENTS': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'ANIMAL_ANESTHESIA': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'ANIMAL_ACCLIMATION_DURATION': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'ANIMAL_FASTING': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'ANIMAL_ENDP_EUTHANASIA': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'ANIMAL_ENDP_TISSUE_COLL_LIST': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'ANIMAL_ENDP_TISSUE_PROC_METHOD': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'ANIMAL_ENDP_CLINICAL_SIGNS': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'HUMAN_FASTING': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'HUMAN_ENDP_CLINICAL_SIGNS': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'CELL_STORAGE': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'CELL_GROWTH_CONTAINER': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'CELL_GROWTH_CONFIG': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'CELL_GROWTH_RATE': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'CELL_INOC_PROC': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'CELL_MEDIA': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'CELL_ENVIR_COND': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'CELL_HARVESTING': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'PLANT_GROWTH_SUPPORT': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'PLANT_GROWTH_LOCATION': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'PLANT_PLOT_DESIGN': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'PLANT_LIGHT_PERIOD': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'PLANT_HUMIDITY': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'PLANT_TEMP': {'type': 'string', **_create_unit_regex_and_message(['°C', 'C'], True)},
                'PLANT_WATERING_REGIME': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'PLANT_NUTRITIONAL_REGIME': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'PLANT_ESTAB_DATE': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'PLANT_HARVEST_DATE': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'PLANT_GROWTH_STAGE': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'PLANT_METAB_QUENCH_METHOD': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'PLANT_HARVEST_METHOD': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'PLANT_STORAGE': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'CELL_PCT_CONFLUENCE': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'CELL_MEDIA_LASTCHANGED': {'type': 'string', 'not':{'enum': NA_VALUES}}},
 'required': ['TREATMENT_SUMMARY'],
 'additionalProperties': False}

sampleprep_schema = \
{'type': 'object',
 'properties': {'SAMPLEPREP_SUMMARY': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'SAMPLEPREP_PROTOCOL_ID': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'SAMPLEPREP_PROTOCOL_FILENAME': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'SAMPLEPREP_PROTOCOL_COMMENTS': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'PROCESSING_METHOD': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'PROCESSING_STORAGE_CONDITIONS': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'EXTRACTION_METHOD': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'EXTRACT_CONCENTRATION_DILUTION': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'EXTRACT_ENRICHMENT': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'EXTRACT_CLEANUP': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'EXTRACT_STORAGE': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'SAMPLE_RESUSPENSION': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'SAMPLE_DERIVATIZATION': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'SAMPLE_SPIKING': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'ORGAN': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'ORGAN_SPECIFICATION': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'CELL_TYPE': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'SUBCELLULAR_LOCATION': {'type': 'string', 'not':{'enum': NA_VALUES}}},
 'required': ['SAMPLEPREP_SUMMARY'],
 'additionalProperties': False}

chromatography_schema = \
{'type': 'object',
 'properties': {'CHROMATOGRAPHY_SUMMARY': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'CHROMATOGRAPHY_TYPE': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'INSTRUMENT_NAME': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'COLUMN_NAME': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'FLOW_GRADIENT': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'FLOW_RATE': {'type': 'string', **_create_unit_regex_and_message(['mL/min', 'uL/min', 'μL/min'], True)},
                'COLUMN_TEMPERATURE': {'type': 'string', **_create_unit_regex_and_message(['°C', 'C'], True)},
                'METHODS_FILENAME': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'SAMPLE_INJECTION': {'type': 'string', **_create_unit_regex_and_message(['μL', 'uL'])},
                'SOLVENT_A': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'SOLVENT_B': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'METHODS_ID': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'COLUMN_PRESSURE': {'type': 'string', **_create_unit_regex_and_message(['psi', 'bar'], True)},
                'INJECTION_TEMPERATURE': {'type': 'string', **_create_unit_regex_and_message(['°C', 'C'], True)},
                'INTERNAL_STANDARD': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'INTERNAL_STANDARD_MT': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'RETENTION_INDEX': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'RETENTION_TIME': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'SAMPLING_CONE': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'ANALYTICAL_TIME': {'type': 'string', **_create_unit_regex_and_message(['min'], True)},
                'CAPILLARY_VOLTAGE': {'type': 'string', **_create_unit_regex_and_message(['V', 'kV'])},
                'MIGRATION_TIME': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'OVEN_TEMPERATURE': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'PRECONDITIONING': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'RUNNING_BUFFER': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'RUNNING_VOLTAGE': {'type': 'string', **_create_unit_regex_and_message(['V', 'kV'])},
                'SHEATH_LIQUID': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'TIME_PROGRAM': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'TRANSFERLINE_TEMPERATURE': {'type': 'string', **_create_unit_regex_and_message(['°C', 'C'])},
                'WASHING_BUFFER': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'WEAK_WASH_SOLVENT_NAME': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'WEAK_WASH_VOLUME': {'type': 'string', **_create_unit_regex_and_message(['μL', 'uL'])},
                'STRONG_WASH_SOLVENT_NAME': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'STRONG_WASH_VOLUME': {'type': 'string', **_create_unit_regex_and_message(['μL', 'uL'])},
                'TARGET_SAMPLE_TEMPERATURE': {'type': 'string', **_create_unit_regex_and_message(['°C', 'C'])},
                'SAMPLE_LOOP_SIZE': {'type': 'string', **_create_unit_regex_and_message(['μL', 'uL'])},
                'SAMPLE_SYRINGE_SIZE': {'type': 'string', **_create_unit_regex_and_message(['μL', 'uL'])},
                'RANDOMIZATION_ORDER': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'CHROMATOGRAPHY_COMMENTS': {'type': 'string', 'not':{'enum': NA_VALUES}}},
 'required': ['CHROMATOGRAPHY_TYPE',
              'INSTRUMENT_NAME',
              'COLUMN_NAME',
              'FLOW_GRADIENT',
              'FLOW_RATE',
              'COLUMN_TEMPERATURE',
              'SOLVENT_A',
              'SOLVENT_B'],
 'additionalProperties': False}

# Add in special case for 'room temperature'.
inj_pattern = chromatography_schema['properties']['INJECTION_TEMPERATURE']['pattern']
new_pattern = '^(' + inj_pattern[1:-1] + ')|(?i:room temperature)$'
chromatography_schema['properties']['INJECTION_TEMPERATURE']['pattern'] = new_pattern
new_message = (' should be a number or range (ex. "5-6") followed by a space with a unit '
               '(ex. "5 V") from the following list: [\'°C\', \'C\'] or "room temperature". '
               'Ignore this when more complicated descriptions are required.')
chromatography_schema['properties']['INJECTION_TEMPERATURE']['pattern_custom_message'] = new_message

analysis_schema = \
{'type': 'object',
 'properties': {'ANALYSIS_TYPE': {'type': 'string', 'enum': ['MS', 'NMR']},
                'LABORATORY_NAME': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'ACQUISITION_DATE': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'SOFTWARE_VERSION': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'OPERATOR_NAME': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'DETECTOR_TYPE': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'ANALYSIS_PROTOCOL_FILE': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'ACQUISITION_PARAMETERS_FILE': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'PROCESSING_PARAMETERS_FILE': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'DATA_FORMAT': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'ACQUISITION_ID': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'ACQUISITION_TIME': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'ANALYSIS_COMMENTS': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'ANALYSIS_DISPLAY': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'INSTRUMENT_NAME': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'INSTRUMENT_PARAMETERS_FILE': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'NUM_FACTORS': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'NUM_METABOLITES': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'PROCESSED_FILE': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'RANDOMIZATION_ORDER': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'RAW_FILE': {'type': 'string', 'not':{'enum': NA_VALUES}}},
 'required': ['ANALYSIS_TYPE'],
 'additionalProperties': False}

results_file_schema = \
{'type': 'object',
 'properties': {'filename': {'type': 'string', 'not':{'enum': NA_VALUES}},
  'UNITS': {'type': 'string', 'not':{'enum': NA_VALUES}},
  'Has m/z': {'type': 'string', 'not':{'enum': NA_VALUES}},
  'Has RT': {'type': 'string', 'not':{'enum': NA_VALUES}},
  'RT units': {'type': 'string', 'not':{'enum': NA_VALUES}}},
 'required': ['filename', 'UNITS'],
 'additionalProperties': False}

ms_schema = \
{'type': 'object',
 'properties': {'INSTRUMENT_NAME': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'INSTRUMENT_TYPE': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'MS_TYPE': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'ION_MODE': {'type': 'string', 'pattern': r'(?i)^(positive|negative|positive, negative|unspecified)$',
                           'pattern_custom_message': (' should be one of "Positive", "Negative", "Positive, Negative", or "Unspecified". '
                                                      'Ignore this when more complicated descriptions are required.')},
                'CAPILLARY_TEMPERATURE': {'type': 'string', **_create_unit_regex_and_message(['°C', 'C'], True)},
                'CAPILLARY_VOLTAGE': {'type': 'string', **_create_unit_regex_and_message(['V', 'kV'])},
                'COLLISION_ENERGY': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'COLLISION_GAS': {'type': 'string', 'pattern': r'(?i)^(nitrogen|argon)$',
                                  'pattern_custom_message': (' should be one of "Nitrogen" or "Argon". '
                                                             'Ignore this when more complicated descriptions are required.')},
                'DRY_GAS_FLOW': {'type': 'string', **_create_unit_regex_and_message(['L/hr', 'L/min'])},
                'DRY_GAS_TEMP': {'type': 'string', **_create_unit_regex_and_message(['°C', 'C'])},
                'FRAGMENT_VOLTAGE': {'type': 'string', **_create_unit_regex_and_message(['V'])},
                'FRAGMENTATION_METHOD': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'GAS_PRESSURE': {'type': 'string', **_create_unit_regex_and_message(['psi', 'psig', 'bar', 'kPa'])},
                'HELIUM_FLOW': {'type': 'string', **_create_unit_regex_and_message(['mL/min'])},
                'ION_SOURCE_TEMPERATURE': {'type': 'string', **_create_unit_regex_and_message(['°C', 'C'])},
                'ION_SPRAY_VOLTAGE': {'type': 'string', **_create_unit_regex_and_message(['V', 'kV'])},
                'IONIZATION': {'type': 'string', 
                               'not':{'oneOf':[{'enum': NA_VALUES}, 
                                               {'pattern': '(?i)^(pos|neg|positive|negative|postive|both)$',
                                                'pattern_custom_message': (' should not be "positive" or "negative". '
                                                                           '"ION_MODE" is where that should be indicated.')}]}},
                'IONIZATION_ENERGY': {'type': 'string', **_create_unit_regex_and_message(['eV'])},
                'IONIZATION_POTENTIAL': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'MASS_ACCURACY': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'PRECURSOR_TYPE': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'REAGENT_GAS': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'SOURCE_TEMPERATURE': {'type': 'string', **_create_unit_regex_and_message(['°C', 'C'])},
                'SPRAY_VOLTAGE': {'type': 'string', **_create_unit_regex_and_message(['kV'])},
                'ACTIVATION_PARAMETER': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'ACTIVATION_TIME': {'type': 'string', **_create_unit_regex_and_message(['ms'])},
                'ATOM_GUN_CURRENT': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'AUTOMATIC_GAIN_CONTROL': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'BOMBARDMENT': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'CDL_SIDE_OCTOPOLES_BIAS_VOLTAGE': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'CDL_TEMPERATURE': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'DATAFORMAT': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'DESOLVATION_GAS_FLOW': {'type': 'string', **_create_unit_regex_and_message(['L/hr', 'L/min'])},
                'DESOLVATION_TEMPERATURE': {'type': 'string', **_create_unit_regex_and_message(['°C', 'C'])},
                'INTERFACE_VOLTAGE': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'IT_SIDE_OCTOPOLES_BIAS_VOLTAGE': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'LASER': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'MATRIX': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'NEBULIZER': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'OCTPOLE_VOLTAGE': {'type': 'string', **_create_unit_regex_and_message(['V'])},
                'PROBE_TIP': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'RESOLUTION_SETTING': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'SAMPLE_DRIPPING': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'SCAN_RANGE_MOVERZ': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'SCANNING': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'SCANNING_CYCLE': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'SCANNING_RANGE': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'SKIMMER_VOLTAGE': {'type': 'string', **_create_unit_regex_and_message(['V'])},
                'TUBE_LENS_VOLTAGE': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'MS_COMMENTS': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'MS_RESULTS_FILE': results_file_schema},
 'required': ['INSTRUMENT_NAME', 'INSTRUMENT_TYPE', 'MS_TYPE', 'ION_MODE'],
 'additionalProperties': False}

nmr_schema = \
{'type': 'object',
 'properties': {'INSTRUMENT_NAME': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'INSTRUMENT_TYPE': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'NMR_EXPERIMENT_TYPE': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'NMR_COMMENTS': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'FIELD_FREQUENCY_LOCK': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'STANDARD_CONCENTRATION': {'type': 'string', **_create_unit_regex_and_message(['mM'])},
                'SPECTROMETER_FREQUENCY': {'type': 'string', **_create_unit_regex_and_message(['MHz'])},
                'NMR_PROBE': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'NMR_SOLVENT': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'NMR_TUBE_SIZE': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'SHIMMING_METHOD': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'PULSE_SEQUENCE': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'WATER_SUPPRESSION': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'PULSE_WIDTH': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'POWER_LEVEL': {'type': 'string', **_create_unit_regex_and_message(['W', 'dB'])},
                'RECEIVER_GAIN': {'type': 'string', **_create_num_regex_and_message(False)},
                'OFFSET_FREQUENCY': {'type': 'string', **_create_unit_regex_and_message(['ppm', 'Hz'])},
                'PRESATURATION_POWER_LEVEL': {'type': 'string', **_create_unit_regex_and_message(['W', 'dB'])},
                'CHEMICAL_SHIFT_REF_CPD': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'TEMPERATURE': {'type': 'string', **_create_unit_regex_and_message(['°C', 'C', 'K'])},
                'NUMBER_OF_SCANS': {'type': 'string', **_create_int_regex_and_message(False)},
                'DUMMY_SCANS': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'ACQUISITION_TIME': {'type': 'string', **_create_unit_regex_and_message(['s'])},
                'RELAXATION_DELAY': {'type': 'string', **_create_unit_regex_and_message(['s', 'ms', 'us', 'μs'])},
                'SPECTRAL_WIDTH': {'type': 'string', **_create_unit_regex_and_message(['ppm', 'Hz'])},
                'NUM_DATA_POINTS_ACQUIRED': {'type': 'string', **_create_int_regex_and_message(False)},
                'REAL_DATA_POINTS': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'LINE_BROADENING': {'type': 'string', **_create_unit_regex_and_message(['Hz'])},
                'ZERO_FILLING': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'APODIZATION': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'BASELINE_CORRECTION_METHOD': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'CHEMICAL_SHIFT_REF_STD': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'BINNED_INCREMENT': {'type': 'string', **_create_unit_regex_and_message(['ppm'])},
                'BINNED_DATA_NORMALIZATION_METHOD': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'BINNED_DATA_PROTOCOL_FILE': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'BINNED_DATA_CHEMICAL_SHIFT_RANGE': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'BINNED_DATA_EXCLUDED_RANGE': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'NMR_RESULTS_FILE': results_file_schema},
 'required': ['INSTRUMENT_NAME',
              'INSTRUMENT_TYPE',
              'NMR_EXPERIMENT_TYPE',
              'SPECTROMETER_FREQUENCY'],
 'additionalProperties': False}

# The 'properties' and 'required' keys are commented out because they are validated in separate 
# functions. Doing it here causes many more spurious errors to be printed.
data_schema = \
{'type': 'array',
 'items': {'type': 'object',
           # 'properties': {'Metabolite': {'type': 'string', 'not':{'enum': METABOLITE_NA_VALUES}}},
           # 'required': ['Metabolite'],
           'additionalProperties': True}}

extended_schema = \
{'type': 'array',
 'items': {'type': 'object',
           # 'properties': {'Metabolite': {'type': 'string', 'not':{'enum': METABOLITE_NA_VALUES}},
           #                'sample_id': {'type': 'string', 'not':{'enum': NA_VALUES}}},
           # 'required': ['Metabolite', 'sample_id'],
           'additionalProperties': True}}

binned_data_schema = \
{'type': 'array',
 'items': {'type': 'object',
           # 'properties': {'Bin range(ppm)': {'type': 'string', 'not':{'enum': NA_VALUES}}},
           # 'required': ['Bin range(ppm)'],
           'additionalProperties': True}}

ms_metabolite_data_schema = \
{'type': 'object',
 'properties': {'Units': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'Data': data_schema,
                'Metabolites': data_schema,
                'Extended': extended_schema},
 'required': ['Units', 'Data'],
 'additionalProperties': False}

nmr_binned_data_schema = \
{'type': 'object',
 'properties': {'Units': {'type': 'string', 'not':{'enum': NA_VALUES}},
                'Data': binned_data_schema},
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
ms_required_schema['then'] = {'properties': {'MS':{'required':['MS_RESULTS_FILE'], 
                                                   'required_message': ('Error: There must be either a "MS_METABOLITE_DATA" '
                                                                        'section or a "MS_RESULTS_FILE" subsection in the '
                                                                        '"MS" section. Neither were found.')}}}

nmr_required_schema = deepcopy(base_required_schema)
nmr_required_schema['properties']['NM'] = nmr_schema
nmr_required_schema['properties']['NMR_METABOLITE_DATA'] = ms_metabolite_data_schema
nmr_required_schema['properties']['NMR_BINNED_DATA'] = nmr_binned_data_schema
nmr_required_schema['required'].extend(['NM'])
nmr_required_schema['if'] = {'allOf':[{'properties': {'NMR_METABOLITE_DATA':{'not':{}}}}, {'properties': {'NMR_BINNED_DATA':{'not':{}}}}]}
nmr_required_schema['then'] = {'properties': {'NM':{'required':['NMR_RESULTS_FILE'],
                                                    'required_message': ('Error: There must be either a "NMR_METABOLITE_DATA" '
                                                                         'section, a "NMR_BINNED_DATA" section or a '
                                                                         '"NMR_RESULTS_FILE" subsection in the '
                                                                         '"NM" section. Neither were found.')}}}

# oneOf prevents specific error messages from being printed, so instead 
# determine whether MS or NMR should be validated and validate using the appropriate schema.
# compiled_schema = \
# {'type': 'object',
#  'oneOf':[ms_required_schema,
#           nmr_required_schema]}
    
    
