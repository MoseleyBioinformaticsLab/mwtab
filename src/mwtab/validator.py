#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
mwtab.validator
~~~~~~~~~~~~~~~

This module contains routines to validate consistency of the ``mwTab``
formatted files, e.g. make sure that ``Samples`` and ``Factors``
identifiers are consistent across the file, and make sure that all
required key-value pairs are present.
"""

from datetime import datetime
from re import match
import io
import sys
import traceback
from collections.abc import Iterable

import jsonschema
import pandas

from .mwschema import ms_required_schema, nmr_required_schema, METABOLITE_NA_VALUES
from .mwschema import NA_VALUES as SCHEMA_NA_VALUES

import mwtab
from mwtab import metadata_column_matching

column_finders = metadata_column_matching.column_finders
NA_VALUES = metadata_column_matching.NA_VALUES

VERBOSE = False
LOG = None

VALIDATION_LOG_HEADER = \
"""Validation Log
{}
mwtab Python Library Version: {}
Source:        {}
Study ID:      {}
Analysis ID:   {}
File format:   {}"""


SUFFIXES = {1: 'st', 2: 'nd', 3: 'rd'}
def ordinal_suffix(num: int) -> str:
    """Return the ordinal suffix for the given integer.
    """
    if 10 <= num % 100 <= 20:
        suffix = 'th'
    else:
        suffix = SUFFIXES.get(num % 10, 'th')
    return str(num) + suffix

def format_column_name(column_name: str, position: int) -> str:
    """Build better message for error printing.
    
    When a column name has duplicates it is more complicated to identify in an 
    error message, so we build a more complex description. Duplicate names should 
    end in a pattern like {{{_\\d_}}}. We can pull out the number and identify 
    which duplicate column (enumerating from the left) is the one needing a 
    description.
    
    Args:
        column_name: The name of the column to build a message for.
        position: A number indicating the position of the column when enumerating from the left most column.
                  Should start from 1.
    
    Returns:
        A string message describing where the column is.
    """
    if column_name.endswith('}}}'):
        new_name = match(mwtab.duplicates_dict.DUPLICATE_KEY_REGEX, column_name).group(1)
        duplicate_num = match(r'.*\{\{\{_(\d+)_\}\}\}', column_name).group(1)
        duplicate_suffix = ordinal_suffix(int(duplicate_num)+1)
        message = f'The {duplicate_suffix} "{new_name}" column at position {position}'
    else:
        message = f'The "{column_name}" column at position {position}'
    
    return message
        


def validate_sub_section_uniqueness(mwtabfile):
    """Validate that the sub-sections don't repeat.
    
    :param mwtabfile: Instance of :class:`~mwtab.mwtab.MWTabFile`.
    :type mwtabfile: :class:`~mwtab.mwtab.MWTabFile`
    """
    errors = []
    if mwtabfile._duplicate_sub_sections:
        for section_name, sub_section in mwtabfile._duplicate_sub_sections.items():
            for sub_section_name in sub_section:
                errors.append({'message': "Error: The section, " + section_name + ", has a sub-section, "
                              + sub_section_name + ", that is duplicated.",
                              'tags': ['consistency'], 'section': section_name, 'sub-section': sub_section_name,
                              'ID': '1', 'name': 'Duplicate Sub-section'})
    return errors


def validate_header_lengths(mwtabfile):
    """Validate that the headers are the correct size.
    
    :param mwtabfile: Instance of :class:`~mwtab.mwtab.MWTabFile`.
    :type mwtabfile: :class:`~mwtab.mwtab.MWTabFile`
    """
    errors = []
    if mwtabfile._short_headers:
        for section in mwtabfile._short_headers:
            errors.append({'message': "Error: The section, " + section + ", has a mismatch between the "
                          "number of headers and the number of elements in each "
                          "line. Either a line(s) has more values than headers or "
                          "there are too few headers.",
                          'tags': ['consistency'], 'section': section, 'sub-section': None,
                          'ID': '2', 'name': 'Bad Headers'})
    return errors


def validate_factors(mwtabfile):
    """Validate that the subject sample factors and factors in the metabolites data match.
    
    :param mwtabfile: Instance of :class:`~mwtab.mwtab.MWTabFile`.
    :type mwtabfile: :class:`~mwtab.mwtab.MWTabFile`
    """
    errors = []
    if mwtabfile._factors:
        factors_dict_1 = {sample: factors for sample, factors in mwtabfile._factors.items()}
        factors_dict_2 = {i["Sample ID"]: i["Factors"] for i in mwtabfile["SUBJECT_SAMPLE_FACTORS"] if i["Sample ID"] in factors_dict_1}
        if factors_dict_1 != factors_dict_2:
            errors.append({'message': "Error: The factors in the METABOLITE_DATA section "
                          "and SUBJECT_SAMPLE_FACTORS section do not match.",
                          'tags': ['consistency'], 'section': 'SUBJECT_SAMPLE_FACTORS', 'sub-section': None,
                          'ID': '3', 'name': 'Factor Mismatch'})
    return errors


def validate_subject_samples_factors(mwtabfile):
    """Validate ``SUBJECT_SAMPLE_FACTORS`` section.

    :param mwtabfile: Instance of :class:`~mwtab.mwtab.MWTabFile`.
    :type mwtabfile: :class:`~mwtab.mwtab.MWTabFile`
    """
    if mwtabfile._input_format == 'mwtab':
        location = 'SUBJECT_SAMPLE_FACTORS entry #{}'
        index_corrector = lambda x: x+1
    else:
        location = 'The SSF at ["SUBJECT_SAMPLE_FACTORS"][{}]'
        index_corrector = lambda x: x
        
    subject_samples_factors_errors = list()
    seen_samples = set()
    for index, subject_sample_factor in enumerate(mwtabfile["SUBJECT_SAMPLE_FACTORS"]):
        if subject_sample_factor["Sample ID"]:
            if subject_sample_factor["Sample ID"] in seen_samples:
                subject_samples_factors_errors.append(
                    {'message': f"Warning: {location.format(index_corrector(index))} has a duplicate Sample ID.",
                     'tags': ['value'], 'section': 'SUBJECT_SAMPLE_FACTORS', 'sub-section': None, 'ID': '4', 'name': 'Duplicate Sample ID in SSF'}
                )
            seen_samples.add(subject_sample_factor["Sample ID"])
        if subject_sample_factor.get("Factors"):           
            duplicate_keys = [re_match.group(1) for key in subject_sample_factor["Factors"]
                              if (re_match := match(mwtab.duplicates_dict.DUPLICATE_KEY_REGEX, key))]
        
            if duplicate_keys:
                subject_samples_factors_errors.append(
                    {'message': f"Warning: {location.format(index_corrector(index))} has the "
                    "following duplicate keys in its Factors:\n\t" + 
                    "\n\t".join(f'"{value}"' for value in duplicate_keys),
                    'tags': ['value'], 'section': 'SUBJECT_SAMPLE_FACTORS', 'sub-section': None, 'ID': '5', 'name': 'Duplicate Factors in SSF'})
        
        if subject_sample_factor.get("Additional sample data"):
            duplicate_keys = [re_match.group(1) for key in subject_sample_factor["Additional sample data"]
                              if (re_match := match(mwtab.duplicates_dict.DUPLICATE_KEY_REGEX, key))]
        
            if duplicate_keys:
                subject_samples_factors_errors.append(
                    {'message': f"Warning: {location.format(index_corrector(index))} has the "
                    "following duplicate keys in its Additional sample data:\n\t" + 
                    "\n\t".join(f'"{value}"' for value in duplicate_keys),
                    'tags': ['value'], 'section': 'SUBJECT_SAMPLE_FACTORS', 'sub-section': None, 'ID': '6', 'name': 'Duplicate Additional Data'})
    return subject_samples_factors_errors


def validate_data(mwtabfile, data_section_key, mwtabfile_tables):
    """Validates ``MS_METABOLITE_DATA``, ``NMR_METABOLITE_DATA``, and ``NMR_BINNED_DATA`` sections.

    :param mwtabfile: Instance of :class:`~mwtab.mwtab.MWTabFile`.
    :type mwtabfile: :class:`~mwtab.mwtab.MWTabFile`
    :param data_section_key: Section key (either MS_METABOLITE_DATA, NMR_METABOLITE_DATA, or NMR_BINNED_DATA)
    :type data_section_key: :py:class:`str`
    :param mwtabfile_tables: Dictionary where the keys are table names and the values are the tables as pandas DataFrames.
    :type mwtabfile_tables: :py:class:`dict`
    """
    data_errors = list()

    subject_sample_factors_sample_id_set = {subject_sample_factor["Sample ID"] for subject_sample_factor in mwtabfile["SUBJECT_SAMPLE_FACTORS"]}
    if mwtabfile._samples:
        data_sample_id_set = set(mwtabfile._samples)
    elif mwtabfile[data_section_key].get("Data"):
        data_keys = list(mwtabfile[data_section_key]["Data"][0].keys())
        data_sample_id_set = set(data_keys[1:])
    else:
        data_sample_id_set = set()

    if mwtabfile._input_format == 'mwtab':
        location = data_section_key
    else:
        location = f'["{data_section_key}"]["Data"]'
    
    if data_sample_id_set - subject_sample_factors_sample_id_set:
        data_errors.append({'message': 'Error: SUBJECT_SAMPLE_FACTORS section missing sample ID(s). '
                           'The following IDs were found in the {} section but not in the SUBJECT_SAMPLE_FACTORS:\n\t{}'.format(
            location,
            "\n\t".join(f'"{value}"' for value in sorted(list(data_sample_id_set - subject_sample_factors_sample_id_set))),
            ),
            'tags': ['consistency'], 'section': 'SUBJECT_SAMPLE_FACTORS', 'sub-section': None, 'ID': '7', 'name': 'Missing Sample ID(s) in SSF'})
    
    # Check if there are duplicate sample names.
    if mwtabfile._samples and (len(mwtabfile._samples) > len(set(mwtabfile._samples))):
        data_errors.append({'message': "Warning: There are duplicate samples in the "
                           f"{location} section.",
                           'tags': ['value'], 'section': data_section_key, 'sub-section': 'Data',
                           'ID': '8', 'name': 'Duplicate Samples in DATA'})
    
    # Check whether mwtabolites in Data are in the Metabolites section.
    if mwtabfile._input_format == 'mwtab':
        metabolites_location = 'METABOLITES'
        data_location = data_section_key
    else:
        metabolites_location = f'["{data_section_key}"]["Metabolites"]'
        data_location = f'["{data_section_key}"]["Data"]'
    
    # If either table does not have a 'Metabolite' column, then a message will be 
    # printed about not being able to find any metabolite. There is a separate 
    # message about not having the column doen elsewhere, so don't run the 
    # 'Metabolite in DATA not METABOLITES' check unless both tables have then column.
    both_have_metabolites_column = True
    metabolites_in_data_section = pandas.Series()
    if mwtabfile_tables['Data'] is not None and not mwtabfile_tables['Data'].empty:
        if 'Metabolite' in mwtabfile_tables['Data'].columns:
            metabolites_in_data_section = mwtabfile_tables['Data'].loc[:, 'Metabolite'].str.strip()
        else:
            both_have_metabolites_column = False
    
    metabolites_in_metabolites = pandas.Series()
    if mwtabfile_tables['Metabolites'] is not None and not mwtabfile_tables['Metabolites'].empty:
        if 'Metabolite' in mwtabfile_tables['Metabolites'].columns:
            metabolites_in_metabolites = mwtabfile_tables['Metabolites'].loc[:, 'Metabolite'].str.strip()
        else:
            both_have_metabolites_column = False
    
    if both_have_metabolites_column:
        data_not_in_met_mask = ~metabolites_in_data_section.isin(metabolites_in_metabolites)
        data_not_in_met = metabolites_in_data_section[data_not_in_met_mask]
        if len(data_not_in_met) > 0 and 'BINNED' not in data_section_key:
            message = (f'Error: The following metabolites in the, {data_location} table '
                      f'were not found in the {metabolites_location} table:\n\t')
            message = message + '\n\t'.join(f'"{value}"' for value in data_not_in_met.values)
            data_errors.append({'message': message, 'tags': ['consistency'], 'section': data_section_key, 'sub-section': 'Data',
                                'ID': '9', 'name': 'Metabolite(s) in DATA not METABOLITES'})
    
    if metabolites_in_data_section.isna().any() or metabolites_in_data_section.isin(METABOLITE_NA_VALUES).any():
        message = f'Error: A metabolite without a name was found in the {data_location} table.'
        data_errors.append({'message': message, 'tags': ['value'], 'section': data_section_key, 'sub-section': 'Data',
                            'ID': '10', 'name': 'Blank Metabolite(s) in DATA'})
    
    duplicate_metabolites_mask = metabolites_in_data_section.duplicated()
    if duplicate_metabolites_mask.any():
        duplicate_metabolites = metabolites_in_data_section[duplicate_metabolites_mask]
        message = f'Warning: The following metabolites in the {data_location} table appear more than once in the table:\n\t'
        message = message + '\n\t'.join(f'"{value}"' for value in duplicate_metabolites.values)
        data_errors.append({'message': message, 'tags': ['value'], 'section': data_section_key, 'sub-section': 'Data',
                            'ID': '11', 'name': 'Duplicate Metabolite(s) in DATA'})

    return data_errors


def validate_metabolites(mwtabfile, data_section_key, mwtabfile_tables):
    """Validate ``METABOLITES`` section.

    :param mwtabfile: Instance of :class:`~mwtab.mwtab.MWTabFile`.
    :type mwtabfile: :class:`~mwtab.mwtab.MWTabFile`
    :param data_section_key: Section key (either MS_METABOLITE_DATA, NMR_METABOLITE_DATA, or NMR_BINNED_DATA)
    :type data_section_key: :py:class:`str`
    :param mwtabfile_tables: Dictionary where the keys are table names and the values are the tables as pandas DataFrames.
    :type mwtabfile_tables: :py:class:`dict`
    """
    implied_pairs = metadata_column_matching.implied_pairs
    
    metabolites_errors = list()
    
    if mwtabfile._input_format == 'mwtab':
        metabolites_location = 'METABOLITES'
        data_location = data_section_key
    else:
        metabolites_location = f'["{data_section_key}"]["Metabolites"]'
        data_location = f'["{data_section_key}"]["Data"]'
    
    # If either table does not have a 'Metabolite' column, then a message will be 
    # printed about not being able to find any metabolite. There is a separate 
    # message about not having the column doen elsewhere, so don't run the 
    # 'Metabolite in METABOLITES not DATA' check unless both tables have then column.
    both_have_metabolites_column = True
    metabolites_in_data_section = pandas.Series()
    if mwtabfile_tables['Data'] is not None and not mwtabfile_tables['Data'].empty:
        if 'Metabolite' in mwtabfile_tables['Data'].columns:
            metabolites_in_data_section = mwtabfile_tables['Data'].loc[:, 'Metabolite'].str.strip()
        else:
            both_have_metabolites_column = False
    
    metabolites_in_metabolites = pandas.Series()
    if mwtabfile_tables['Metabolites'] is not None and not mwtabfile_tables['Metabolites'].empty:
        if 'Metabolite' in mwtabfile_tables['Metabolites'].columns:
            metabolites_in_metabolites = mwtabfile_tables['Metabolites'].loc[:, 'Metabolite'].str.strip()
        else:
            both_have_metabolites_column = False
        
    
    if both_have_metabolites_column:
        met_not_in_data_mask = ~metabolites_in_metabolites.isin(metabolites_in_data_section)
        met_not_in_data = metabolites_in_metabolites[met_not_in_data_mask]
        if len(met_not_in_data) > 0:
            message = (f'Error: The following metabolites in the {metabolites_location} table '
                      f'were not found in the {data_location} table:\n\t')
            message = message + '\n\t'.join(f'"{value}"' for value in met_not_in_data.values)
            metabolites_errors.append({'message': message, 'tags': ['consistency'], 'section': data_section_key, 'sub-section': 'Metabolites',
                                       'ID': '12', 'name': 'Metabolite(s) in METABOLITES not DATA'})
    
    if metabolites_in_metabolites.isna().any() or metabolites_in_metabolites.isin(METABOLITE_NA_VALUES).any():
        message = f'Error: A metabolite without a name was found in the {metabolites_location} table.'
        metabolites_errors.append({'message': message, 'tags': ['value'], 'section': data_section_key, 'sub-section': 'Metabolites',
                                   'ID': '13', 'name': 'Blank Metabolite(s) in METABOLITES'})
    
    duplicate_metabolites_mask = metabolites_in_metabolites.duplicated()
    if duplicate_metabolites_mask.any():
        duplicate_metabolites = metabolites_in_metabolites[duplicate_metabolites_mask]
        message = f'Warning: The following metabolites in the {metabolites_location} table appear more than once in the table:\n\t'
        message = message + '\n\t'.join(f'"{value}"' for value in duplicate_metabolites.values)
        metabolites_errors.append({'message': message, 'tags': ['value'], 'section': data_section_key, 'sub-section': 'Metabolites',
                                   'ID': '14', 'name': 'Duplicate Metabolite(s) in METABOLITES'})

    # Check if fields/columns are recognized variations and report the standardized name to the user.
    df = mwtabfile_tables['Metabolites']
    columns = {column:column.lower().strip() for column in df.columns}
    found_columns = {}
    columns_to_standard_columns = {}
    for name, finder in column_finders.items():
        if column_matches := finder.name_dict_match(columns):
            found_columns[name] = column_matches
            if name not in df.columns:
                for column_name in column_matches:
                    # We ignore if the column name is only off due to capitalization because lowering is common and easy.
                    if column_name.lower() != name:
                        message = (f'Warning: {format_column_name(column_name, df.columns.get_loc(column_name)+1)} '
                                   f'in the {metabolites_location} table, '
                                   f'matches a standard column name, "{name}". '
                                   'If this match was not in error, the column should be renamed to '
                                   'the standard name or a name that doesn\'t resemble the standard name.')
                        metabolites_errors.append({'message': message, 'tags': ['value'], 'section': data_section_key, 'sub-section': 'Metabolites',
                                                   'ID': '15', 'name': 'Standard Column Name Match'})
            for column_name in column_matches:
                if column_name in columns_to_standard_columns:
                    columns_to_standard_columns[column_name].append(name)
                else:
                    columns_to_standard_columns[column_name] = [name]
                
                value_mask = finder.values_series_match(df.loc[:, column_name].astype('string[pyarrow]'), na_values = NA_VALUES)
                if not value_mask.all():
                    message = (f'Warning: {format_column_name(column_name, df.columns.get_loc(column_name)+1)} '
                               f'in the {metabolites_location} table, '
                               f'matches a standard column name, "{name}", '
                               'and some of the values in the column do not match the expected type or format for that column. '
                               'The non-matching values are:\n')
                    message += df.loc[~value_mask, column_name].to_string()
                    metabolites_errors.append({'message': message, 'tags': ['value'], 'section': data_section_key, 'sub-section': 'Metabolites',
                                               'ID': '16', 'name': 'METABOLITES Bad Standard Values'})
    
    # When certain columns are found in METABOLITES, look for the implied pair and warn if it isn't there. 
    # For example, other_id and other_id_type and retention_index and retention_index_type.
    # Also check that they both have data in the same rows.
    for name, matches in found_columns.items():
        if name in implied_pairs:
            implied_pairs_not_found = [column for column in implied_pairs[name] if column not in found_columns]
            for column in implied_pairs_not_found:
                message = (f'Warning: The column "{matches[0]}" was found in the {metabolites_location} table, '
                           f'but this column implies that another column, "{column}", '
                           'should also exist, and that column was not found.')
                metabolites_errors.append({'message': message, 'tags': [], 'section': data_section_key, 'sub-section': 'Metabolites',
                                           'ID': '17', 'name': 'Missing Implied Column'})
            
            implied_pairs_found = [column for column in implied_pairs[name] if column in found_columns]
            for column in implied_pairs_found:
                parent_mask = (df.loc[:, matches[0]].isna() | df.loc[:, matches[0]].isin(NA_VALUES))
                child_mask = (df.loc[:, found_columns[column][0]].isna() | df.loc[:, found_columns[column][0]].isin(NA_VALUES))
                if (parent_mask != child_mask).any():
                    message = (f'Warning: The column pair, "{matches[0]}" and "{found_columns[column][0]}", '
                               'in the METABOLITES table should have data in the '
                               'same rows, but at least one row has data in one '
                               'column and nothing in the other.')
                    metabolites_errors.append({'message': message, 'tags': [], 'section': data_section_key, 'sub-section': 'Metabolites', 
                                               'ID': '18', 'name': 'Paired Columns Value Mismatch'})
    
    # If the other_id column is found, print message about making individual database ID columns.
    # I thought about checking to see if there were database IDs in the column first, but I'm not sure if it's worth the effort.
    if 'other_id' in found_columns:
        message = ('Warning: The standard column, "other_id", was '
                   f'found in the METABOLITES table as "{found_columns["other_id"][0]}". '
                   'If this column contains database IDs for standard databases such '
                   'as KEGG, PubChem, HMDB, etc., it is recommended to make individual '
                   'columns for these and not lump them together into a less descriptive '
                   '"other_id" column.')
        metabolites_errors.append({'message': message, 'tags': [], 'section': data_section_key, 'sub-section': 'Metabolites', 
                                   'ID': '19', 'name': '"other_id" Column'})
    
    # If a column in df matches multiple standard names, print a warning to prefer separating them.
    for column_name, standard_names in columns_to_standard_columns.items():
        if len(standard_names) > 1:
            message = (f'Warning: The column, "{column_name}", in the {metabolites_location} table '
                       f'was matched to multiple standard names, {standard_names}. This is a good indication '
                       'that the values in that column should be split into the appropriate '
                       'individual columns.')
            metabolites_errors.append({'message': message, 'tags': [], 'section': data_section_key, 'sub-section': 'Metabolites',
                                       'ID': '20', 'name': 'Multiple Standard Name Match'})

    return metabolites_errors


def validate_extended(mwtabfile, data_section_key, mwtabfile_tables):
    """Validate ``EXTENDED_MS_METABOLITE_DATA``, ``EXTENDED_NMR_METABOLITE_DATA``, and ``EXTENDED_NMR_BINNED_DATA`` sections.

    :param mwtabfile: Instance of :class:`~mwtab.mwtab.MWTabFile`.
    :type mwtabfile: :class:`~mwtab.mwtab.MWTabFile`
    :param data_section_key: Section key (either MS_METABOLITE_DATA, NMR_METABOLITE_DATA, or NMR_BINNED_DATA)
    :type data_section_key: :py:class:`str`
    :param mwtabfile_tables: Dictionary where the keys are table names and the values are the tables as pandas DataFrames.
    :type mwtabfile_tables: :py:class:`dict`
    """
    extended_errors = list()

    sample_id_set = {subject_sample_factor["Sample ID"] for subject_sample_factor in
                     mwtabfile["SUBJECT_SAMPLE_FACTORS"]}
    
    if mwtabfile._input_format == 'mwtab':
        extended_location = 'EXTENDED_METABOLITE_DATA'
        ssf_string = 'in the SUBJECT_SAMPLE_FACTORS section'
        
    else:
        extended_location = f' in ["{data_section_key}"]["Extended"]'
        ssf_string = 'in ["SUBJECT_SAMPLE_FACTORS"]'
    
    df = mwtabfile_tables['Extended']
    if "sample_id" not in df.columns:
        message = f"Error: The {extended_location} table does not have a column for \"sample_id\"."
        extended_errors.append({'message': message, 'tags': ['format'], 'section': data_section_key, 'sub-section': 'Extended',
                                'ID': '21', 'name': 'Missing "sample_id" in EXTENDED'})
    else:
        extended_id_set = set(df.loc[:, 'sample_id'])
        not_in_ssf = extended_id_set - sample_id_set
        if not_in_ssf:
            message = (f"Error: The {extended_location} table has Sample IDs that were not found "
                       f"{ssf_string}. Those IDs are:\n\t" + '\n\t'.join(f'"{value}"' for value in not_in_ssf))
            extended_errors.append({'message': message, 'tags': ['consistency'], 'section': data_section_key, 'sub-section': 'Extended',
                                    'ID': '22', 'name': 'Missing Sample ID(s) in EXTENDED'})
        
        if df.loc[:, 'sample_id'].isna().any() or df.loc[:, 'sample_id'].isin(SCHEMA_NA_VALUES).any():
            message = f'Error: A Sample ID without a name was found in the {extended_location} table.'
            extended_errors.append({'message': message, 'tags': ['value'], 'section': data_section_key, 'sub-section': 'Extended',
                                    'ID': '35', 'name': 'Blank Sample ID(s) in EXTENDED'})
    
    if 'Metabolite' in df.columns and (df.loc[:, 'Metabolite'].isna().any() or df.loc[:, 'Metabolite'].isin(METABOLITE_NA_VALUES).any()):
        message = f'Error: A metabolite without a name was found in the {extended_location} table.'
        extended_errors.append({'message': message, 'tags': ['value'], 'section': data_section_key, 'sub-section': 'Extended',
                                   'ID': '36', 'name': 'Blank Metabolite(s) in EXTENDED'})

    return extended_errors


def validate_metabolite_names(mwtabfile, data_section_key):
    """Validate that a header didn't get identified as a metabolite.

    :param mwtabfile: Instance of :class:`~mwtab.mwtab.MWTabFile`.
    :type mwtabfile: :class:`~mwtab.mwtab.MWTabFile`
    :param data_section_key: Section key (either MS_METABOLITE_DATA, NMR_METABOLITE_DATA, or NMR_BINNED_DATA)
    :type data_section_key: :py:class:`str`
    """
    # Lowered versions of headers to look for, and common mispellings.
    list_headers = ["samples", "factors", "bin range(ppm)", "metabolite_name", "metabolite name"]
    
    if "Metabolites" in mwtabfile[data_section_key]:
        metabolites_section_bad_names = [data_dict["Metabolite"] for data_dict in mwtabfile[data_section_key]["Metabolites"] 
                                         if 'Metabolite' in data_dict and data_dict["Metabolite"].lower() in list_headers]
    else:
        metabolites_section_bad_names = []
    
    if "Data" in mwtabfile[data_section_key]:
        data_section_bad_names = [data_dict["Metabolite"] for data_dict in mwtabfile[data_section_key]["Data"] 
                                  if 'Metabolite' in data_dict and data_dict["Metabolite"].lower() in list_headers]
    else:
        data_section_bad_names = []
    
    if "Extended" in mwtabfile[data_section_key]:
        extended_section_bad_names = [data_dict["Metabolite"] for data_dict in mwtabfile[data_section_key]["Extended"] 
                                      if 'Metabolite' in data_dict and data_dict["Metabolite"].lower() in list_headers]
    else:
        extended_section_bad_names = []
    
    if mwtabfile._input_format == 'mwtab':
        desc_string = 'in the {} table'
        metabolite_string = 'METABOLITES'
        data_string = 'METABOLITE_DATA' if not 'BINNED' in data_section_key else 'BINNED_DATA'
        extended_string = 'EXTENDED_METABOLITE_DATA'
    else:
        desc_string = f'in the table at ["{data_section_key}"]{{}}'
        metabolite_string = '[\"Metabolites\"]'
        data_string = '[\"Data\"]'
        extended_string = '[\"Extended\"]'
    errors = []
    message = ("Warning: There is a metabolite name, "
               f"\"{{}}\", {desc_string} that is probably wrong. "
               "It is close to a header name and is likely due to a badly constructed Tab file.")
    for name in metabolites_section_bad_names:
        errors.append({'message': message.format(name, metabolite_string), 
                       'tags': ['value'], 'section': data_section_key, 'sub-section': 'Metabolites',
                       'ID': '23', 'name': 'Bad Metabolite Name'})
    for name in data_section_bad_names:
        errors.append({'message': message.format(name, data_string), 
                       'tags': ['value'], 'section': data_section_key, 'sub-section': 'Data',
                       'ID': '23', 'name': 'Bad Metabolite Name'})
    for name in extended_section_bad_names:
        errors.append({'message': message.format(name, extended_string), 
                       'tags': ['value'], 'section': data_section_key, 'sub-section': 'Extended',
                       'ID': '23', 'name': 'Bad Metabolite Name'})
    
    return errors
    


def _format_section_name(section, subsection):
    section_name = section
    if 'METABOLITE_DATA' in section:
        if subsection == 'Extended':
            section_name = 'EXTENDED_' + section
        elif subsection == 'Metabolites':
            section_name = 'METABOLITES'
    return section_name
        

def _gen_jsonschema_format_string(error: jsonschema.ValidationError, 
                                  error_path_len: int, 
                                  schema: dict, 
                                  mwtabfile: 'mwtab.mwtab.MWTabFile') -> [str, bool|None, str|None, str|None]:
    """Generates a few variables that create_better_error_messages will use in creating better error messages.
    
    Args:
        error: The error to create messages for.
        error_path_len: The len of error.relative_path.
        schema: The jsonschema given to create_better_error_messages.
        mwtabfile: The mwtabfile that was validated against.
    
    Returns:
        A tuple where the first element is a format string used in messages, 
        the second is a boolean that is True if the key the error is about is required,
        the third element is the section of the mwTab file the error came from 
        (None when this does not make sense), and 
        the fourth element is the subsection of the mwTab file the error came from 
        (None when this does not make sense).
        The boolean will be None when key requirement does not make any sense.
    """
    format_string = ''
    key_is_required = None
    section = None
    subsection = None
    # This should be an additional properties error.
    if error_path_len == 0:
        pass
    # 4 levels deep should be an issue with a metabolite in MS_METABOLITE_DATA or something in SUBJECT_SAMPLE_FACTORS.
    elif error_path_len == 4:
        section = error.relative_path[0]
        subsection = error.relative_path[1]
        element = error.relative_path[2]
        result_key = error.relative_path[3]
                
        if isinstance(subsection, int):
            if mwtabfile._input_format == 'mwtab':
                format_string = f'for the "{result_key}" in "{element}" in entry {subsection+1} of the "{_format_section_name(section, subsection)}" section'
            else:
                format_string = f'in ["{section}"][{subsection}]["{element}"]["{result_key}"]'
            
            # This is likely a key in the Factors or Additional data dicts of SUBJECT_SAMPLE_FACTORS, so no keys are required.
            required = []
        else:
            if mwtabfile._input_format == 'mwtab':
                format_string = f'in the "{_format_section_name(section, subsection)}" section, in row number {element+1}, for the "{result_key}" column'
            else:
                format_string = f'in ["{section}"]["{subsection}"][{element}]["{result_key}"]'
            
            required = schema['properties'][section]['properties'][subsection]['items'].get('required')
        key_is_required = required and result_key in required
        
    # MS_RESULTS_FILE or NMR_RESULTS_FILE or SUBJECT_SAMPLE_FACTORS.
    elif error_path_len == 3:
        section = error.relative_path[0]
        subsection = error.relative_path[1]
        result_key = error.relative_path[2]
        
        if isinstance(subsection, int):
            if mwtabfile._input_format == 'mwtab':
                format_string = f'for the "{result_key}" in entry {subsection+1} of the "{_format_section_name(section, subsection)}" section'
            else:
                format_string = f'in ["{section}"][{subsection}]["{result_key}"]'
            
            required = schema['properties'][section]['items'].get('required')
        else:
            if mwtabfile._input_format == 'mwtab':
                format_string = f'for the subsection, "{subsection}", in the "{_format_section_name(section, subsection)}" section, for the "{result_key}" attribute'
            else:
                format_string = f'in ["{section}"]["{subsection}"]["{result_key}"]'
            
            required = schema['properties'][section]['properties'][subsection].get('required')
        key_is_required = required and result_key in required
        
    elif error_path_len == 1:
        section = error.relative_path[0]
        
        if mwtabfile._input_format == 'mwtab':
            format_string = f'in the "{section}" section '
        else:
            format_string = f'in ["{section}"]'
    else:
        section = error.relative_path[0]
        subsection = error.relative_path[1]
        
        if mwtabfile._input_format == 'mwtab':
            if section == 'METABOLOMICS WORKBENCH':
                format_string = f'for "{subsection}" in the file header'
            else:
                format_string = f'for the subsection, "{subsection}", in the "{_format_section_name(section, subsection)}" section'
        else:
            format_string = f'in ["{section}"]["{subsection}"]'
        
        required = schema['properties'][section].get('required')
        key_is_required = required and subsection in required
        
    return format_string, key_is_required, section, subsection
    
VALIDATOR_TO_TAG = {
    'additionalProperties': ['format'],
    'minProperties': ['value'],
    'required': ['format'],
    'minLength': ['value'],
    'maxLength': ['value'],
    'minItems': ['value'],
    'type': ['format'],
    'enum': ['value'],
    'format': ['value'],
    'pattern': ['value'],
    'minimum': ['value'],
    'maximum': ['value'],
    'uniqueItems': ['value'],
    'not': ['value']
    }
def create_better_error_messages(errors_generator: Iterable[jsonschema.exceptions.ValidationError], 
                                 mwtabfile: 'mwtab.mwtab.MWTabFile',
                                 schema: dict) -> list[str]:
    """Create better error messages for jsonschema validation errors.
        
    Args:
        errors_generator: The generator returned from validator.iter_errors().
        mwtabfile: The mwtabfile that was validated against.
        schema: The schema used to create errors_generator. Used to get more information and craft better messages.
    
    Returns:
        A list of the errors encountered.
    """
    
    # try:
    #     jsonschema.validate(instance=instance, schema=schema)
    # except jsonschema.ValidationError as e:
    #     ## code to easily see the contents of the error for building a better message.
    #     for key, value in e._contents().items():
    #         print(key, value)
    #         print()
    
    # for key, value in error._contents().items():
    #     print(key, value)
    #     print()
    
    errors = []
    for error in errors_generator:
        validator = error.validator
        validator_value = error.validator_value
        if prefix := error.schema.get(f'{validator}_prefix'):
            message = prefix
        else:
            message = "Error: "
        custom_message = ""
        
        
        error_path_len = len(error.relative_path)
        format_string, key_is_required, section, subsection = _gen_jsonschema_format_string(error, error_path_len, schema, mwtabfile)
        
        # There is one special case where 'not' is not only being used for {'not': {'enum':NA_VALUES}}.
        # This is for IONIZATION in the MS section. This code is a more generalized approach to determine 
        # which schema within a 'oneOf' is the one that actually triggered the validation error. It 
        # is probably overkill for this special case, but might be relevant in the future.
        if validator == "not" and 'oneOf' in error.schema['not']:
            real_schema_found = False
            i = 0
            while not real_schema_found and i < len(error.schema['not']['oneOf']):
                try:
                    jsonschema.validate(error.instance, {'not':error.schema['not']['oneOf'][i]})
                except jsonschema.ValidationError as error2:
                    real_schema_found = True
                    custom_message_keys = [key for key in error2.schema['not'] if 'custom_message' in key]
                    error.schema = error2.schema['not']
                    if custom_message_keys:
                        validator = match(r'(.*)_custom_message', custom_message_keys[0]).group(1)
                    else:
                        # Assuming that the schema is only a single keyword.
                        validator = list(error2.schema['not'].keys())[0]
                i += 1
        
        
        if custom_message_attr := error.schema.get(f'{validator}_custom_message'):
            custom_message = custom_message_attr
        elif message_attr := error.schema.get(f'{validator}_message'):
            message = message_attr
        elif validator == "additionalProperties":
            if 'were' in error.message:
                bad_keys = match(r'.* \((.*) were unexpected\).*', error.message).group(1)
                bad_keys = bad_keys.replace("'", '"')
                if error_path_len == 0:
                    message = message + f'Unknown or invalid sections, {bad_keys}.'
                else:
                    message = message + f'Unknown or invalid subsections, {bad_keys}, {format_string}.'
            else:
                bad_key = match(r'.* \(\'(.*)\' was unexpected\).*', error.message).group(1)
                if error_path_len == 0:
                    message = message + f'Unknown or invalid section, "{bad_key}".'
                else:
                    message = message + f'Unknown or invalid subsection, "{bad_key}", {format_string}.'
        elif validator == "minProperties":
            custom_message = " cannot be empty."
        elif validator == "required":
            required_property = match(r"\'(.*)\'", error.message).group(1)
            message += f'The required property, "{required_property}", {format_string} is missing.'
        # Commenting these out for now because they aren't used in mwschema and I don't want to test them.
        # elif validator == "dependencies":
        #     message += "The entry " + "[%s]" % "][".join(repr(index) for index in error.relative_path) + " is missing a dependent property.\n"
        #     message += error.message
        # elif validator == "dependentRequired":
        #     message += "The entry " + "[%s]" % "][".join(repr(index) for index in error.relative_path) + " is missing a dependent property.\n"
        #     message += error.message
        elif validator == "minLength":
            if validator_value == 1 and isinstance(error.instance, str):
                custom_message = " cannot be an empty string."
            else:
                custom_message = " is too short."
        elif validator == "maxLength":
            custom_message = " is too long."
        elif validator == "minItems":
            if validator_value == 1:
                custom_message = " cannot be empty."
            else:
                custom_message = " must have at least " + str(validator_value) + " items."
        elif validator == "type":
            if type(validator_value) == list:
                custom_message = " is not any of the allowed types: ["
                for allowed_type in validator_value:
                    custom_message += "\'" + allowed_type + "\', "
                custom_message = custom_message[:-2]
                custom_message += "]."
            else:
                custom_message = " is not of type \"" + validator_value + "\"."
        elif validator == "enum":
            custom_message = " is not one of [" + "%s" % ", ".join(repr(index) for index in validator_value) + "]."
        elif validator == "format":
            custom_message = " is not a valid " + validator_value + "."
        elif validator == "pattern":
            custom_message = " does not match the regular expression pattern " + str(validator_value)
        elif validator == "minimum":
            custom_message = " must be greater than or equal to " + str(validator_value) + "."
        elif validator == "maximum":
            custom_message = " must be less than or equal to " + str(validator_value) + "."
        elif validator == "uniqueItems":
            custom_message = " has non-unique elements."
        # For now 'not' is assumed to only be for {'not': {'enum':NA_VALUES}}.
        # If a new use case for 'not' gets added, see if adding 'not_custom_message' 
        # like what is done for Factors in SUBJECT_SAMPLE_FACTORS will do what is needed.
        elif validator == 'not':
            message = message + f'An empty value or a null value was detected {format_string}.'
            
            if mwtabfile._input_format == 'mwtab':
                container_noun = 'subsection'
            else:
                container_noun = 'key'
            
            if key_is_required:
                message = message + f' A legitimate value should be provided for this required {container_noun}.'
            else:
                message = message + f' Either a legitimate value should be provided for this {container_noun}, or it should be removed altogether.'
        else:
            message += error.message
        
        
        if custom_message:
            str_instance = str(error.instance)
            if len(str_instance) < 50:
                message = message + f"The value, \"{str_instance}\", {format_string}" + custom_message
            else:
                message = message + f"The value {format_string}" + custom_message
        errors.append({'message': message, 'tags': VALIDATOR_TO_TAG.get(validator, []), 
                       'section': section, 'sub-section': subsection,
                       'ID': '24', 'name': 'JSON Schema Error: ' + validator})
    return errors


def validate_schema(mwtabfile, schema):
    """Validate section of ``mwTab`` formatted file.
    
    :param mwtabfile: Instance of :class:`~mwtab.mwtab.MWTabFile`.
    :type mwtabfile: :class:`~mwtab.mwtab.MWTabFile`
    :param schema: JSON Schema schema.
    :type schema: :py:class:`dict`
    :return: JSON Schema errors.
    :rtype: :py:class:`list`
    """
    validator = jsonschema.validators.validator_for(schema)
    format_checker = jsonschema.FormatChecker()
    validator = validator(schema=schema, format_checker=format_checker)
    # errors_generator = validator.iter_errors(mwtabfile)
    return create_better_error_messages(validator.iter_errors(mwtabfile), mwtabfile, schema)
    
  
    
def validate_table_values(mwtabfile, data_section_key, mwtabfile_tables, na_values = NA_VALUES):
    """Validate the values of all table sections.

    :param mwtabfile: Instance of :class:`~mwtab.mwtab.MWTabFile`.
    :type mwtabfile: :class:`~mwtab.mwtab.MWTabFile`
    :param data_section_key: Section key (either MS_METABOLITE_DATA, NMR_METABOLITE_DATA, or NMR_BINNED_DATA)
    :type data_section_key: :py:class:`str`
    :param mwtabfile_tables: Dictionary where the keys are table names and the values are the tables as pandas DataFrames.
    :type mwtabfile_tables: :py:class:`dict`
    :param na_values: List of values to consider null values. 
    :type na_values: :py:class:`list`
    """
    if mwtabfile._input_format == 'mwtab':
        data_location = data_section_key
        met_location = 'METABOLITES'
        extended_location = 'EXTENDED_METABOLITE_DATA'
        
        header_names = {
            'Data': 'Bin range(ppm)' if 'BINNED' in data_section_key else 'Samples',
            'Metabolites': 'metabolite_name',
            'Extended': 'metabolite_name'
            }
    else:
        data_location = f'["{data_section_key}"]["Data"]'
        met_location = f'["{data_section_key}"]["Metabolites"]'
        extended_location = f'["{data_section_key}"]["Extended"]'
        
        header_names = {
            'Data': 'Bin range(ppm)' if 'BINNED' in data_section_key else 'Metabolite',
            'Metabolites': 'Metabolite',
            'Extended': 'Metabolite'
            }
    
    message_strings = {
        'Data': data_location,
        'Metabolites': met_location,
        'Extended': extended_location
        }
    headers = {
        'Data': mwtabfile._raw_binned_header if 'BINNED' in data_section_key else mwtabfile._raw_samples,
        'Metabolites': mwtabfile._raw_metabolite_header,
        'Extended': mwtabfile._raw_extended_metabolite_header
        }
    errors = []
    for table_name in mwtabfile.table_names:
        if table_name in mwtabfile[data_section_key]:
            if headers[table_name] and header_names[table_name] not in headers[table_name]:
                message = (f"Error: The {message_strings[table_name]} table does not have a column for \"{header_names[table_name]}\". "
                           "It is likely misspelled or using a common incorrect substitute.")
                errors.append({'message': message, 'tags': ['format'], 'section': data_section_key, 'sub-section': table_name,
                               'ID': '34', 'name': 'Missing Header'})
                        
            records = mwtabfile[data_section_key][table_name]
            if len(records) > 0:
                columns = list(records[0].keys())
                if not all([list(data_dict.keys()) == columns for data_dict in records]):
                    message = (f'Error: The {message_strings[table_name]} table '
                              'does not have the same columns for every row.')
                    errors.append({'message': message, 'tags': ['consistency'], 'section': data_section_key, 'sub-section': table_name,
                                   'ID': '25', 'name': 'Inconsistent Columns'})
            
            data_df = mwtabfile_tables[table_name]
            
            # Look for empty column names.
            if any(name == '' for name in data_df.columns):
                message = (f'Error: Column(s) with no name were found in the {message_strings[table_name]} table.')
                errors.append({'message': message, 'tags': ['value'], 'section': data_section_key, 'sub-section': table_name,
                               'ID': '26', 'name': 'Column With No Name'})
            
            # Look for NA values that aren't the empty string.
            # Hunter did not like this in validation since normalizing NA values is not too difficult as a cleaning step.
            # na_values_sans_empty_string = [value for value in na_values if value != '']
            # values_mask = data_df.isin(na_values_sans_empty_string)
            # if values_mask.any().any():
            #     na_values_in_df = [f'"{value}"' for value in pandas.unique(data_df[values_mask].values.ravel()) if not pandas.isna(value)]
            #     message = (f'Warning: Nonstandard NA values were found in the {message_strings[table_name]} table. '
            #                'NA values should be expressed using the empty string unless '
            #                'the empty string means something else, such as "not found", '
            #                'in which case there should be both NA and NF values and no empty string. '
            #                'The following values were the nonstandard values found:\n\t')
            #     message = message + '\n\t'.join(na_values_in_df)
            #     errors.append(message)
            
            # Look for completely null columns.
            null_columns = data_df.isna().all() | data_df.isin(na_values).all()
            null_columns = null_columns[null_columns]
            if len(null_columns) > 0:
                for column in null_columns.index:
                    message = (f'Warning: {format_column_name(column, data_df.columns.get_loc(column)+1)} '
                               f'in the {message_strings[table_name]} table has all null values.')
                    errors.append({'message': message, 'tags': ['value'], 'section': data_section_key, 'sub-section': table_name,
                                   'ID': '27', 'name': 'Null Column'})
            
            # Look for overbalanced values, so if 90% of a column is dominated by a single value print a warning.
            for i, column in enumerate([column for column in data_df.columns if column != 'Metabolite']):
                temp_column = data_df.loc[:, column].astype(str)
                value_counts = temp_column.value_counts(dropna=False)
                value_counts = value_counts / value_counts.sum()
                if any(value_counts > .9) and len(value_counts) > 1:
                    message = (f'Warning: {format_column_name(column, data_df.columns.get_loc(column)+1)} '
                               f'in the {message_strings[table_name]} table may have incorrect values. '
                               '90% or more of the values are the same, but 10% or less are different.')
                    errors.append({'message': message, 'tags': ['value'], 'section': data_section_key, 'sub-section': table_name,
                                   'ID': '28', 'name': 'Possible Bad Column Values'})
            
            # Look for duplicate rows.
            if data_df.duplicated().any():
                message = f"Warning: There are duplicate rows in the {message_strings[table_name]} table."
                errors.append({'message': message, 'tags': ['value'], 'section': data_section_key, 'sub-section': table_name,
                               'ID': '29', 'name': 'Duplicate Rows'})
            
            # Look for duplicate column names. We skip Data because there is a separate check to look for duplicate samples.
            if not table_name == 'Data':
                columns = [column if not column.endswith('}}}') else match(mwtab.duplicates_dict.DUPLICATE_KEY_REGEX, column).group(1) 
                                           for column in data_df.columns]
                if len(columns) > len(set(columns)):
                    message = f"Warning: There are duplicate column names in the {message_strings[table_name]} table."
                    errors.append({'message': message, 'tags': ['value'], 'section': data_section_key, 'sub-section': table_name,
                                   'ID': '30', 'name': 'Duplicate Column Names'})
    return errors


def validate_polarity(mwtabfile, data_section_key, mwtabfile_tables):
    """Validate that polarity columns are mono valued.
    
    :param mwtabfile_tables: Dictionary where the keys are table names and the values are the tables as pandas DataFrames.
    :type mwtabfile_tables: :py:class:`dict`
    :param data_section_key: Section key (either MS_METABOLITE_DATA, NMR_METABOLITE_DATA, or NMR_BINNED_DATA)
    :type data_section_key: :py:class:`str`
    :param mwtabfile_tables: Dictionary where the keys are table names and the values are the tables as pandas DataFrames.
    :type mwtabfile_tables: :py:class:`dict`
    """
    if mwtabfile._input_format == 'mwtab':
        location = 'METABOLITES'
    else:
        location = f'["{data_section_key}"]["Metabolites"]'
    
    errors = []    
    df = mwtabfile_tables['Metabolites']
    column_finder = column_finders['polarity']
    columns = {column:column.lower().strip() for column in df.columns}
    if column_matches := column_finder.name_dict_match(columns):
        for column_match in column_matches:
            pos_values = df.loc[:, column_match].str.lower().isin(['pos', 'positive', '+'])
            neg_values = df.loc[:, column_match].str.lower().isin(['neg', 'negative', '+'])
            if pos_values.any() and neg_values.any():
                message = (f'Error: The "{column_match}" column in the {location} table '
                           'indicates multiple polarities in a single analysis, and '
                           'this should not be. A single mwTab file is supposed to be '
                           'restricted to a single analysis. This means multiple MS '
                           'runs under different settings should each be in their own file.')
                errors.append({'message': message, 'tags': ['format'], 'section': data_section_key, 'sub-section': 'Metabolites',
                               'ID': '31', 'name': 'Multiple Polarities'})
                break
    return errors


def validate_file(mwtabfile: 'mwtab.mwtab.MWTabFile', 
                  ms_schema: dict = ms_required_schema,
                  nmr_schema: dict = nmr_required_schema,
                  verbose: bool = False) -> (str, list[dict]):
    """Validate ``mwTab`` formatted file.
    
    Note that some of the validations are pretty strict to account for the majority of cases, 
    but if warranted could be ignored. For example, COLUMN_PRESSURE in the CHROMATOGRAPHY 
    section will print a warning if the value is not a single number or range of numbers 
    followed by a unit, but there might be some situations where the method is complex 
    and thus the column pressure is not static. So something like 
    "60 bar at starting conditions. 180 bar at %A" would be required to accurately 
    describe the COLUMN_PRESSURE, and would be valid. So in these kinds of situations 
    the warning printed can safely be ignored.
    
    Args:
        mwtabfile: The file to be validated.
        ms_schema: jsonschema to validate both the base parts of the file and the MS specific parts of the file.
        nmr_schema: jsonschema to validate both the base parts of the file and the NMR specific parts of the file.
        verbose: whether to be verbose or not.
    
    Returns:
        Error messages as a single string and error messages in JSON form. If verbose is True, then the single string will be None.
    """
    # setup
    if not verbose:
        error_stout = io.StringIO()
    else:
        error_stout = sys.stdout

    # generate validation log header(s)
    file_format = mwtabfile.source.split("/")[-1] if "https://www.metabolomicsworkbench.org/" in mwtabfile.source else \
        mwtabfile.source.split(".")[1]
    print(VALIDATION_LOG_HEADER.format(
        str(datetime.now()),
        mwtab.__version__,
        mwtabfile.source,
        mwtabfile.study_id,
        mwtabfile.analysis_id,
        file_format
    ), file=error_stout)

    # create list to collect validation errors
    errors = list()
    
    # Get tables as dataframes.
    mwtabfile_tables = {}
    for table_name in mwtabfile.table_names:
        mwtabfile_tables[table_name] = mwtabfile.get_table_as_pandas(table_name)
    
    if 'NM' in mwtabfile:
        errors.extend(validate_schema(mwtabfile, nmr_schema))
    else:
        if 'MS' not in mwtabfile:
            message = ('Error: No "MS" or "NM" section was found, '
                       'so analysis type could not be determined. '
                       'Mass spec will be assumed.')
            errors.append({'message': message, 'tags': ['format'], 'section': None, 
                           'sub-section': None, 'ID': '32', 'name': 'No MS or NM Section'})
        errors.extend(validate_schema(mwtabfile, ms_schema))
    
    # validate SUBJECT_SAMPLE_FACTORS
    errors.extend(validate_subject_samples_factors(mwtabfile))
    errors.extend(validate_factors(mwtabfile))

    # validate ..._DATA sections
    data_section_key = mwtabfile.data_section_key
    if data_section_key:
        errors.extend(validate_data(mwtabfile, data_section_key, mwtabfile_tables))

        if data_section_key in ("MS_METABOLITE_DATA", "NMR_METABOLITE_DATA"):
            if "Metabolites" in mwtabfile[data_section_key].keys():
                errors.extend(validate_metabolites(mwtabfile, data_section_key, mwtabfile_tables))
            else:
                if mwtabfile._input_format == 'mwtab':
                    location = 'METABOLITES'
                else:
                    location = f'["{data_section_key}"]["Metabolites"]'
                message = f"Warning: Missing {location} section."
                errors.append({'message': message, 'tags': ['format'], 'section': None, 
                               'sub-section': None, 'ID': '33', 'name': 'Missing METABOLITES Section'})
        
        if "Extended" in mwtabfile[data_section_key].keys():
            errors.extend(validate_extended(mwtabfile, data_section_key, mwtabfile_tables))
        
        errors.extend(validate_metabolite_names(mwtabfile, data_section_key))
        errors.extend(validate_table_values(mwtabfile, data_section_key, mwtabfile_tables))
        errors.extend(validate_polarity(mwtabfile, data_section_key, mwtabfile_tables))
    
    errors.extend(validate_header_lengths(mwtabfile))
    errors.extend(validate_sub_section_uniqueness(mwtabfile))
    

    # finish writing validation/error log
    if errors:
        print("Status: Contains Validation Issues", file=error_stout)
        print("Number of Issues: {}\n".format(len(errors)), file=error_stout)
        tags = []
        warning_count = 0
        error_messages = []
        for error in errors:
            if error['message'].startswith('Warning'):
                warning_count += 1
            tags += error['tags']
            error_messages.append(error['message'])
        tag_counts = pandas.Series(tags).value_counts()
        print("Number of Warnings: {}\n".format(warning_count), file=error_stout)
        print("Number of Value Errors: {}\n".format(tag_counts.get('value', 0)), file=error_stout)
        print("Number of Consistency Errors: {}\n".format(tag_counts.get('consistency', 0)), file=error_stout)
        print("Number of Format Errors: {}\n".format(tag_counts.get('format', 0)), file=error_stout)
        try:
            print("Issue Log:\n" + "\n".join(error_messages), file=error_stout)
        # I have only seen this exception when trying to run from Windows 
        # Command Prompt while also redirecting the output to a file.
        # I have no idea how to get this to trigger through a test in pytest.
        except UnicodeEncodeError:
            try:
                print("An error occurred when trying to print the issue log with unicode. Trying UTF-8.")
                print(("Issue Log:\n" + "\n".join(error_messages)).encode('utf-8'), file=error_stout)
            except Exception as e:
                print("An error occurred when trying to print the issue log, so it could not be printed.")
                print("Issue Log Exception:\n")
                traceback.print_exception(e, file=error_stout)
    else:
        print("Status: Passing", file=error_stout)

    if verbose:
        return None, errors
    else:
        return error_stout.getvalue(), errors


