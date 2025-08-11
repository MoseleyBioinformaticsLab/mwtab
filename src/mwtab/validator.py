#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
mwtab.validator
~~~~~~~~~~~~~~~

This module contains routines to validate consistency of the ``mwTab``
formatted files, e.g. make sure that ``Samples`` and ``Factors``
identifiers are consistent across the file, make sure that all
required key-value pairs are present.
"""

from collections import OrderedDict
from datetime import datetime
from re import match
import io
import sys
import traceback
from collections.abc import Iterable

import jsonschema

from .mwschema import section_schema_mapping, base_schema

import mwtab
from mwtab import metadata_column_matching

column_finders = metadata_column_matching.column_finders

VERBOSE = False
LOG = None

METABOLITES_REGEXS = {
    "hmdb_id": {
        r"(?i)[\s|\S]{,}(HMDB)",
        r"(?i)(Human Metabolome D)[\S]{,}",
    },
    "inchi_key": {
        r"(?i)(inchi)[\S]{,}",
    },
    "kegg_id": {
        r"(?i)(kegg)$",
        r"(?i)(kegg)(\s|_)(i)",
    },
    "moverz": {
        r"(?i)(m/z)",
    },
    "moverz_quant": {
        r"(?i)(moverz)(\s|_)(quant)",
        r"(?i)(quan)[\S]{,}(\s|_)(m)[\S]{,}(z)",
    },
    "other_id": {
        r"(?i)(other)(\s|_)(id)$",
    },
    "other_id_type": {
        r"(?i)(other)(\s|_)(id)(\s|_)(type)$",
    },
    "pubchem_id": {
        r"(?i)(pubchem)[\S]{,}",
    },
    "retention_index": {
        r"(?i)(ri)$",
        r"(?i)(ret)[\s|\S]{,}(index)",
    },
    "retention_index_type": {
        r"(?i)(ri)(\s|_)(type)",
    },
    "retention_time": {
        r"(?i)(r)[\s|\S]{,}(time)[\S]{,}",
    },
}

ITEM_SECTIONS = {
    "METABOLOMICS WORKBENCH",
    "PROJECT",
    "STUDY",
    "ANALYSIS",
    "SUBJECT",
    "COLLECTION",
    "TREATMENT",
    "SAMPLEPREP",
    "CHROMATOGRAPHY",
    "MS",
    "NM",
}

VALIDATION_LOG_HEADER = \
"""Validation Log
{}
mwtab Python Library Version: {}
Source:        {}
Study ID:      {}
Analysis ID:   {}
File format:   {}"""



def validate_sub_section_uniqueness(mwtabfile):
    """Validate that the sub-sections don't repeat.
    
    :param mwtabfile: Instance of :class:`~mwtab.mwtab.MWTabFile`.
    :type mwtabfile: :class:`~mwtab.mwtab.MWTabFile` or
                     :py:class:`collections.OrderedDict`
    """
    errors = []
    if mwtabfile._duplicate_sub_sections:
        for block_name, sub_section in mwtabfile._duplicate_sub_sections.items():
            for sub_section_name in sub_section:
                errors.append("Error: The block, " + block_name + ", has a sub-section, "
                              + sub_section_name + ", that is duplicated.")
    return errors


def validate_header_lengths(mwtabfile):
    """Validate that the headers are the correct size.
    
    :param mwtabfile: Instance of :class:`~mwtab.mwtab.MWTabFile`.
    :type mwtabfile: :class:`~mwtab.mwtab.MWTabFile` or
                     :py:class:`collections.OrderedDict`
    """
    errors = []
    if mwtabfile._short_headers:
        for block in mwtabfile._short_headers:
            errors.append("Error: The block, " + block + ", has a mismatch between the "
                          "number of headers and the number of elements in each "
                          "line. Either a line(s) has more values than headers or "
                          "there are too few headers.")
    return errors


def validate_factors(mwtabfile):
    """Validate that the subject sample factors and factors in the metabolites data match.
    
    :param mwtabfile: Instance of :class:`~mwtab.mwtab.MWTabFile`.
    :type mwtabfile: :class:`~mwtab.mwtab.MWTabFile` or
                     :py:class:`collections.OrderedDict`
    """
    errors = []
    if mwtabfile._factors:
        factors_dict_1 = {sample: factors for sample, factors in mwtabfile._factors.items()}
        factors_dict_2 = {i["Sample ID"]: i["Factors"] for i in mwtabfile["SUBJECT_SAMPLE_FACTORS"] if i["Sample ID"] in factors_dict_1}
        if factors_dict_1 != factors_dict_2:
            errors.append("SUBJECT_SAMPLE_FACTORS: The factors in the "
                          "METABOLITE_DATA section and SUBJECT_SAMPLE_FACTORS section "
                          "do not match.")
    return errors


def validate_metabolite_headers(mwtabfile):
    """Validate that the headers in the METABOLITES, EXTENDED_METABOLITES, and BINNED_DATA section are unique.
    
    :param mwtabfile: Instance of :class:`~mwtab.mwtab.MWTabFile`.
    :type mwtabfile: :class:`~mwtab.mwtab.MWTabFile` or
                     :py:class:`collections.OrderedDict`
    """
    errors = []
    if mwtabfile._metabolite_header and (len(mwtabfile._metabolite_header) > len(set(mwtabfile._metabolite_header))):
        errors.append("METABOLITES: There are duplicate metabolite headers in the "
                      "METABOLITES section. (metabolite_name line)")
    if mwtabfile._extended_metabolite_header and (len(mwtabfile._extended_metabolite_header) > len(set(mwtabfile._extended_metabolite_header))):
        errors.append("EXTENDED_METABOLITES: There are duplicate metabolite headers in the "
                      "EXTENDED_METABOLITES section. (metabolite_name line)")
    if mwtabfile._binned_header and (len(mwtabfile._binned_header) > len(set(mwtabfile._binned_header))):
        errors.append("BINNED_DATA: There are duplicate headers in the "
                      "BINNED_DATA section. (Bin range(ppm) line)")
    return errors


def validate_samples(mwtabfile, data_section_key):
    """Validate that the samples in SUBJECT_SAMPLE_FACTORS and METABOLITE_DATA or BINNED_DATA are the same.
    
    :param mwtabfile: Instance of :class:`~mwtab.mwtab.MWTabFile`.
    :type mwtabfile: :class:`~mwtab.mwtab.MWTabFile` or
                     :py:class:`collections.OrderedDict`
    :param data_section_key: Section key (either MS_METABOLITE_DATA, NMR_METABOLITE_DATA, or NMR_BINNED_DATA)
    :type data_section_key: :py:class:`str`
    """
    errors = []
        
    if mwtabfile._samples and (len(mwtabfile._samples) > len(set(mwtabfile._samples))):
        errors.append("METABOLITE_DATA: There are duplicate Samples in the "
                      "METABOLITE_DATA section. (Samples line)")
    return errors


def validate_subject_samples_factors(mwtabfile):
    """Validate ``SUBJECT_SAMPLE_FACTORS`` section.

    :param mwtabfile: Instance of :class:`~mwtab.mwtab.MWTabFile`.
    :type mwtabfile: :class:`~mwtab.mwtab.MWTabFile` or
                     :py:class:`collections.OrderedDict`
    """
    subject_samples_factors_errors = list()
    
    seen_samples = set()
    for index, subject_sample_factor in enumerate(mwtabfile["SUBJECT_SAMPLE_FACTORS"]):
        if not subject_sample_factor["Subject ID"]:
            subject_samples_factors_errors.append(
                "SUBJECT_SAMPLE_FACTORS: Entry #{} missing Subject ID.".format(index+1)
            )
        if not subject_sample_factor["Sample ID"]:
            subject_samples_factors_errors.append(
                "SUBJECT_SAMPLE_FACTORS: Entry #{} missing Sample ID.".format(index + 1)
            )
        else:
            if subject_sample_factor["Sample ID"] in seen_samples:
                subject_samples_factors_errors.append(
                    "SUBJECT_SAMPLE_FACTORS: Entry #{} has a duplicate Sample ID.".format(index + 1)
                )
            seen_samples.add(subject_sample_factor["Sample ID"])
        if subject_sample_factor.get("Factors"):
            ssf_factors = subject_sample_factor["Factors"]
            for factor_key in ssf_factors:
                if not ssf_factors[factor_key]:
                    subject_samples_factors_errors.append(
                        "SUBJECT_SAMPLE_FACTORS: Entry #{} missing value for Factor {}.".format(index + 1, factor_key)
                    )
            
            duplicate_keys = [match(mwtab.mwtab.DUPLICATE_KEY_REGEX, key).group(1) 
                              for key in subject_sample_factor["Factors"]
                              if match(mwtab.mwtab.DUPLICATE_KEY_REGEX, key)]
        
            if duplicate_keys:
                subject_samples_factors_errors.append("SUBJECT_SAMPLE_FACTORS: Entry #" + str(index + 1) + 
                                                      " has the following duplicate keys in Factors:\n\t" + 
                                                      "\n\t".join(duplicate_keys))
        
        if subject_sample_factor.get("Additional sample data"):
            ssf_asd = subject_sample_factor["Additional sample data"]
            for additional_key in ssf_asd:
                if not ssf_asd[additional_key]:
                    subject_samples_factors_errors.append(
                        "SUBJECT_SAMPLE_FACTORS: Entry #{} missing value for Additional sample data {}.".format(
                            index + 1, additional_key
                        )
                    )
            
            duplicate_keys = [match(mwtab.mwtab.DUPLICATE_KEY_REGEX, key).group(1) 
                              for key in subject_sample_factor["Additional sample data"]
                              if match(mwtab.mwtab.DUPLICATE_KEY_REGEX, key)]
        
            if duplicate_keys:
                subject_samples_factors_errors.append("SUBJECT_SAMPLE_FACTORS: Entry #" + str(index + 1) + 
                                                      " has the following duplicate keys in Additional sample data:\n\t" + 
                                                      "\n\t".join(duplicate_keys))
    return subject_samples_factors_errors


def validate_data(mwtabfile, data_section_key, null_values, metabolites):
    """Validates ``MS_METABOLITE_DATA``, ``NMR_METABOLITE_DATA``, and ``NMR_BINNED_DATA`` sections.

    :param mwtabfile: Instance of :class:`~mwtab.mwtab.MWTabFile`.
    :type mwtabfile: :class:`~mwtab.mwtab.MWTabFile` or
                     :py:class:`collections.OrderedDict`
    :param data_section_key: Section key (either MS_METABOLITE_DATA, NMR_METABOLITE_DATA, or NMR_BINNED_DATA)
    :type data_section_key: :py:class:`str`
    :param bool null_values: whether null values are present.
    :param bool metabolites: whether to check that metabolites in the Data section are in the Metabolites section.
    """
    data_errors = list()

    subject_sample_factors_sample_id_set = {subject_sample_factor["Sample ID"] for subject_sample_factor in mwtabfile["SUBJECT_SAMPLE_FACTORS"]}
    if mwtabfile._samples:
        data_sample_id_set = set(mwtabfile._samples)
    elif mwtabfile[data_section_key]["Data"]:
        data_keys = list(mwtabfile[data_section_key]["Data"][0].keys())
        data_sample_id_set = set(data_keys[1:])
    else:
        data_sample_id_set = set()
    if metabolites and "Metabolites" in mwtabfile[data_section_key]:
        metabolites_in_metabolites_section = [data_dict["Metabolite"] for data_dict in mwtabfile[data_section_key]["Metabolites"]]

    # Removed for mwTab File Spec. 1.5
    # if subject_sample_factors_sample_id_set - data_sample_id_set:
    #     data_errors.append("{}: Section missing data entry for sample(s): {}.".format(
    #         data_section_key,
    #         subject_sample_factors_sample_id_set - data_sample_id_set
    #     ))
    if data_sample_id_set - subject_sample_factors_sample_id_set:
        data_errors.append("SUBJECT_SAMPLE_FACTORS: Section missing sample ID(s). The following IDs were found in the {} section but not in the SUBJECT_SAMPLE_FACTORS: {}".format(
            data_section_key,
            sorted(list(data_sample_id_set - subject_sample_factors_sample_id_set)),
        ))

    for index, metabolite_dict in enumerate(mwtabfile[data_section_key]["Data"]):
        # if set(list(metabolite.keys())[1:]) != subject_sample_factors_sample_id_set:
        #     print(len(subject_sample_factors_sample_id_set), len(metabolite) - 1)
        #     print(
        #         "{}: Metabolite \"{}\" missing data entry for {} samples".format(
        #             data_section_key,
        #             metabolite[list(metabolite.keys())[0]],
        #             len(subject_sample_factors_sample_id_set - set(list(metabolite.keys())[1:]))
        #         ),
        #         file=error_stream
        #     )
        
        # Check whether the metabolite is in the Metabolites section.
        if metabolites and "Metabolites" in mwtabfile[data_section_key]:
            metabolite = metabolite_dict["Metabolite"]
            if metabolite not in metabolites_in_metabolites_section:
                data_errors.append("DATA: Data entry #{}, \"{}\", is not in the Metabolites section.".format(index + 1, metabolite))
        
        # Check if there are null values.
        # I think there are issues with this code. This function is never called with null_values = True, so I'm just going to leave it.
        if null_values:
            for data_point_key in metabolite_dict.keys():
                if data_point_key != "Metabolite":
                    try:
                        float(metabolite_dict[data_point_key])
                    except ValueError:
                        # metabolite_dict[data_point_key] = ""
                        data_errors.append(
                            "{}: Data entry #{} contains non-numeric value converted to \"\".".format(data_section_key, index + 1))

    return data_errors


def validate_metabolites(mwtabfile, data_section_key, mwtabfile_tables):
    """Validate ``METABOLITES`` section.

    :param mwtabfile: Instance of :class:`~mwtab.mwtab.MWTabFile`.
    :type mwtabfile: :class:`~mwtab.mwtab.MWTabFile` or
                     :py:class:`collections.OrderedDict`
    :param data_section_key: Section key (either MS_METABOLITE_DATA, NMR_METABOLITE_DATA, or NMR_BINNED_DATA)
    :type data_section_key: :py:class:`str`
    :param mwtabfile_tables: Dictionary where the keys are table names and the values are the tables as pandas DataFrames.
    :type mwtabfile_tables: :py:class:`dict`
    """
    implied_pairs = metadata_column_matching.implied_pairs
    
    metabolites_errors = list()
    
    metabolites_in_data_section = [data_dict["Metabolite"] for data_dict in mwtabfile[data_section_key]["Data"]]

    for index, metabolite_dict in enumerate(mwtabfile[data_section_key]["Metabolites"]):
        # TODO make this check work from dataframes.
        # Check whether the metabolite is in the Data section.
        metabolite = metabolite_dict["Metabolite"]
        if metabolite not in metabolites_in_data_section:
            metabolites_errors.append("METABOLITES: Data entry #{}, \"{}\", is not in the Data section.".format(index + 1, metabolite))
        
        # TODO remove this check since the one below it should replace it.
        # Check if fields are recognized variations and report the standardized name to the user.
        for field_key in list(metabolite_dict.keys())[1:]:
            if not any(k == field_key for k in METABOLITES_REGEXS.keys()):
                for regex_key in METABOLITES_REGEXS.keys():
                    if any(match(p, field_key) for p in METABOLITES_REGEXS[regex_key]):
                        metabolites_errors.append("METABOLITES: Data entry #{} contains field name \"{}\" which matches a commonly used field name \"{}\".".format(index + 1, field_key, regex_key))
                        field_key = regex_key
                        break
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
                    metabolites_errors.append(f'Warning: The column, "{column_name}", matches a standard column name, "{name}". '
                                              'If this match was not in error, the column should be renamed to '
                                              'the standard name or a name that doesn\'t resemble the standard name.')
            for column_name in column_matches:
                if column_name in columns_to_standard_columns:
                    columns_to_standard_columns[column_name].append(name)
                else:
                    columns_to_standard_columns[column_name] = [name]
                
                value_mask = finder.value_series_match(df.loc[:, column_name].astype('string[pyarrow]'))
                if not value_mask.all():
                    error_message = (f'Warning: The column, "{column_name}", matches a standard column name, "{name}", '
                                    'and some of the values in the column do not match the expected type or format for that column. '
                                    'The non-matching values are:\n')
                    error_message += str(df.loc[~value_mask, column_name])
                    metabolites_errors.append(error_message)
    
    # When certain columns are found in METABOLITES, look for the implied pair and warn if it isn't there. 
    # For example, other_id and other_id_type and retention_index and retention_index_type.
    # Also check that they both have data in the same rows.
    for name, matches in found_columns.items():
        if name in implied_pairs:
            implied_pairs_not_found = [column for column in implied_pairs[name] if column not in found_columns]
            for column in implied_pairs_not_found:
                metabolites_errors.append(f'Warning: The column "{matches[0]}" was found in the METABOLITES table, '
                                          'but this column implies that another column, "{column}", '
                                          'should also exist, and that column was not found.')
            
            implied_pairs_found = [column for column in implied_pairs[name] if column in found_columns]
            for column in implied_pairs_found:
                parent_mask = df.loc[:, name].isna()
                child_mask = df.loc[:, columns].isna()
                if (parent_mask != child_mask).any():
                    metabolites_errors.append(f'Warning: The column pair, "{matches[0]}" and "{found_columns[column][0]}", '
                                              'in the METABOLITES table should have data in the '
                                              'same rows, but at least one row has data in one '
                                              'column and nothing in the other.')
    
    # If the other_id column is found, print message about making individual database ID columns.
    # I thought about checking to see if there were database IDs in the column first, but I'm not sure if it's worth the effort.
    if 'other_id' in found_columns:
        metabolites_errors.append('Warning: The standard column, "other_id", was '
                                  f'found in the METABOLITES table as "{found_columns["other_id"][0]}". '
                                  'If this column contains database IDs for standard databases such '
                                  'as KEGG, PubChem, HMDB, etc., it is recommended to make individual '
                                  'columns for these and not lump them together into a less descriptive '
                                  '"other_id" column.')
    
    # If a column in df matches multiple standard names, print a warning to prefer separating them.
    for column_name, standard_names in columns_to_standard_columns.items():
        if len(standard_names) > 1:
            metabolites_errors.append(f'Warning: The column, "{column_name}", in the METABOLTIES table '
                                      'was matched to multiple standard names. This is a good indication '
                                      'that the values in that column should be split into the appropriate '
                                      'individual columns.')

    return metabolites_errors


def validate_extended(mwtabfile, data_section_key, mwtabfile_tables):
    """Validate ``EXTENDED_MS_METABOLITE_DATA``, ``EXTENDED_NMR_METABOLITE_DATA``, and ``EXTENDED_NMR_BINNED_DATA`` sections.

    :param mwtabfile: Instance of :class:`~mwtab.mwtab.MWTabFile`.
    :type mwtabfile: :class:`~mwtab.mwtab.MWTabFile` or
                     :py:class:`collections.OrderedDict`
    :param data_section_key: Section key (either MS_METABOLITE_DATA, NMR_METABOLITE_DATA, or NMR_BINNED_DATA)
    :type data_section_key: :py:class:`str`
    :param mwtabfile_tables: Dictionary where the keys are table names and the values are the tables as pandas DataFrames.
    :type mwtabfile_tables: :py:class:`dict`
    """
    extended_errors = list()

    sample_id_set = {subject_sample_factor["Sample ID"] for subject_sample_factor in
                     mwtabfile["SUBJECT_SAMPLE_FACTORS"]}
    
    df = mwtabfile_tables['Extended']
    if "sample_id" not in df.columns:
        extended_errors.append("Error: The EXTENDED data table does not have a column for \"sample_id\".")
    else:
        extended_id_set = set(df.loc[:, 'sample_id'])
        not_in_ssf = sample_id_set - extended_id_set
        if not_in_ssf:
            extended_errors.append("Error: The EXTENDED data table has Sample IDs that were not found in the "
                                   "SUBJECT_SAMPLE_FACTORS section. Those IDs are:\n" + '\n'.join(not_in_ssf))

    # TODO confirm the new way identifies the same missing sample IDs.
    # for index, extended_data in enumerate(mwtabfile[data_section_key]["Extended"]):
    #     if "sample_id" not in extended_data.keys():
    #         extended_errors.append("EXTENDED_{}: Data entry #{} missing Sample ID.".format(data_section_key, index + 1))
    #     elif not extended_data["sample_id"] in sample_id_set:
    #         extended_errors.append(
    #             "EXTENDED_{}: Data entry #{} contains Sample ID \"{}\" not found in SUBJECT_SAMPLE_FACTORS section.".format(
    #                 data_section_key, index + 1, extended_data["sample_id"]
    #             ))

    return extended_errors


def validate_metabolite_names(mwtabfile, data_section_key):
    """Validate that a header didn't get identified as a metabolite.

    :param mwtabfile: Instance of :class:`~mwtab.mwtab.MWTabFile`.
    :type mwtabfile: :class:`~mwtab.mwtab.MWTabFile` or
                     :py:class:`collections.OrderedDict`
    :param data_section_key: Section key (either MS_METABOLITE_DATA, NMR_METABOLITE_DATA, or NMR_BINNED_DATA)
    :type data_section_key: :py:class:`str`
    """
    # Lowered versions of headers to look for, and common mispellings.
    list_headers = ["samples", "factors", "bin range(ppm)", "metabolite_name", "metabolite name"]
    
    if "Metabolites" in mwtabfile[data_section_key]:
        metabolites_section_bad_names = [data_dict["Metabolite"] for data_dict in mwtabfile[data_section_key]["Metabolites"] if data_dict["Metabolite"].lower() in list_headers]
    else:
        metabolites_section_bad_names = []
    
    if "Data" in mwtabfile[data_section_key]:
        data_section_bad_names = [data_dict["Metabolite"] for data_dict in mwtabfile[data_section_key]["Data"] if data_dict["Metabolite"].lower() in list_headers]
    else:
        data_section_bad_names = []
    
    if "Extended" in mwtabfile[data_section_key]:
        extended_section_bad_names = [data_dict["Metabolite"] for data_dict in mwtabfile[data_section_key]["Extended"] if data_dict["Metabolite"].lower() in list_headers]
    else:
        extended_section_bad_names = []
    
    errors = []
    message = ("Warning: There is a metabolite name, "
               "\"{}\", in the {} section that is probably wrong. "
               "It is close to a header name and is likely due to a badly constructed Tab file.")
    for name in metabolites_section_bad_names:
        errors.append(message.format(name, "[\"" + data_section_key + "\"][\"Metabolites\"] \ METABOLITES"))
    for name in data_section_bad_names:
        errors.append(message.format(name, "[\"" + data_section_key + "\"][\"Data\"] \ METABOLITE_DATA" if not 'BINNED' in data_section_key else "[\"" + data_section_key + "\"][\"Data\"] \ BINNED_DATA"))
    for name in extended_section_bad_names:
        errors.append(message.format(name, "[\"" + data_section_key + "\"][\"Extended\"] \ EXTENDED_METABOLITE_DATA"))
    
    return errors
    



def print_better_error_messages(errors_generator: Iterable[jsonschema.exceptions.ValidationError], pattern_messages: dict|None = None) -> bool:
    """Print better error messages for jsonschema validation errors.
    
    Args:
        errors_generator: the generator returned from validator.iter_errors().
        pattern_messages: if an error is the pattern type, then look up the section and subsection
          of the mwTab that failed the pattern in this dict and see if there is a custom message.
    
    Returns:
        A list of the errors encountered.
    """
    if not pattern_messages:
        pattern_messages = {}
        
    errors = []
    for error in errors_generator:
        
        message = ""
        custom_message = ""
        
        if error.validator == "minProperties":
            custom_message = " cannot be empty."
        elif error.validator == "required":
            required_property = match(r"(\'.*\')", error.message).group(1)
            if len(error.relative_path) == 0:
                message += "The required property " + required_property + " is missing."
            else:
                message += "The entry " + "[%s]" % "][".join(repr(index) for index in error.relative_path) + " is missing the required property " + required_property + "."
        elif error.validator == "dependencies":
            message += "The entry " + "[%s]" % "][".join(repr(index) for index in error.relative_path) + " is missing a dependent property.\n"
            message += error.message
        elif error.validator == "dependentRequired":
            message += "The entry " + "[%s]" % "][".join(repr(index) for index in error.relative_path) + " is missing a dependent property.\n"
            message += error.message
        elif error.validator == "minLength":
            if error.validator_value == 1 and isinstance(error.instance, str):
                custom_message = " cannot be an empty string."
            else:
                custom_message = " is too short."
        elif error.validator == "maxLength":
            custom_message = " is too long."
        elif error.validator == "minItems":
            if error.validator_value == 1:
                custom_message = " cannot be empty."
            else:
                custom_message = " must have at least " + str(error.validator_value) + " items."
        elif error.validator == "type":
            if type(error.validator_value) == list:
                custom_message = " is not any of the allowed types: ["
                for allowed_type in error.validator_value:
                    custom_message += "\'" + allowed_type + "\', "
                custom_message = custom_message[:-2]
                custom_message += "]."
            else:
                custom_message = " is not of type \"" + error.validator_value + "\"."
        elif error.validator == "enum":
            custom_message = " is not one of [" + "%s" % ", ".join(repr(index) for index in error.validator_value) + "]."
        elif error.validator == "format":
            custom_message = " is not a valid " + error.validator_value + "."
        elif error.validator == "pattern":
            if (section := error.relative_path[-2]) in pattern_messages and (subsection := error.relative_path[-1]) in pattern_messages[section]:
                custom_message = pattern_messages[section][subsection]
            else:
                custom_message = " does not match the regular expression pattern " + str(error.validator_value)
        elif error.validator == "minimum":
            custom_message = " must be greater than or equal to " + str(error.validator_value) + "."
        elif error.validator == "maximum":
            custom_message = " must be less than or equal to " + str(error.validator_value) + "."
        elif error.validator == "uniqueItems":
            custom_message = " has non-unique elements."
        else:
            errors.append(error)
            # print(error, file=sys.stderr)
        
        
        if custom_message:
            message = message + "The value for " + "[%s]" % "][".join(repr(index) for index in error.relative_path) + custom_message
        errors.append(message)
        # print("Error:  " + message, file=sys.stderr)
    return errors


def validate_section_schema(mwtabfile, section, schema, section_key):
    """Validate section of ``mwTab`` formatted file.
    
    :param mwtabfile: Instance of :class:`~mwtab.mwtab.MWTabFile`.
    :type mwtabfile: :class:`~mwtab.mwtab.MWTabFile` or
                     :py:class:`collections.OrderedDict`
    :param section: Section of :class:`~mwtab.mwtab.MWTabFile`.
    :type section: :py:class:`collections.OrderedDict`
    :param schema: Schema definition.
    :type schema: :py:class:`~schema.schema`
    :param str section_key: Section key.
    :return: Validated section.
    :rtype: :py:class:`collections.OrderedDict`
    """
    validator = jsonschema.validators.validator_for(schema)
    format_checker = jsonschema.FormatChecker()
    validator = validator(schema=schema, format_checker=format_checker)
    print_better_error_messages(validator.iter_errors(mwtabfile))
    
    
    
    
    
    mwtabfile = {'PROJECT': {}}
    dict_for_Schema = OrderedDict()
    for section_key, section in mwtabfile.items():
        dict_for_Schema[section_key] = section
    errors = []
    try:
        base_schema.validate(dict_for_Schema)
    except Exception as e:
        messages = []
        for auto in e.autos:
            print(auto)
            if auto:
                if auto.startswith("Wrong keys"):
                    bad_key = match(r"Wrong keys ('.*') in .*", auto).group(1)
                    auto = "Unknown or malformed keys " + bad_key + "."
                elif auto.startswith("Wrong key"):
                    bad_key = match(r"Wrong key ('.*') .*", auto).group(1)
                    auto = "Unknown or malformed key " + bad_key + "."
                messages.append(auto)
        errors.append("SCHEMA: There is an issue with the top level of the data.\n" + 
                      "\n".join(messages))
        # errors.append("SCHEMA: There is an issue with the top level of the data.\n" + 
        #               "\n".join([message for message in e.autos[1:] if message]))
    # Need to manually check this because Schema can't do conditional checks.
    if 'MS' in mwtabfile and 'CHROMATOGRAPHY' not in mwtabfile:
        if errors[-1].startswith("SCHEMA: There is an issue with the top level of the data."):
            errors[-1] = errors[-1] + "\nMissing key: 'CHROMATOGRAPHY'"
        else:
            errors.append("SCHEMA: There is an issue with the top level of the data.\nMissing key: 'CHROMATOGRAPHY'")

    # validate PROJECT, STUDY, ANALYSIS... and Schemas
    # TODO see if this catches multiple errors in a single section. For example, bad SUBJECT_TYPE and bad SUBJECT_SPECIES in the SUBJECT section.
    mwtabfile = {'MS': {'Data':2, 'Units':3}}
    errors = []
    for section_key, section in mwtabfile.items():
        try:
            schema = section_schema_mapping[section_key]
            # section = validate_section_schema(section, schema, section_key, error_stout)
            section, schema_errors = validate_section_schema(section, schema, section_key)
            errors.extend(schema_errors)
            # mwtabfile[section_key] = section
        except Exception as e:
            errors.append("SCHEMA: Section \"{}\" does not match the allowed schema. ".format(section_key) + str(e))
    
    
    
    
    
    
    
    
    
    
    schema_errors = list()

    keys_to_delete = []
    if section_key in ITEM_SECTIONS:
        for key in section.keys():
            if not section[key]:
                schema_errors.append("{}: Contains item \"{}\" with null value.".format(section_key, key))
                keys_to_delete.append(key)

    return schema.validate(section), schema_errors
    
def validate_table_values(mwtabfile, data_section_key, mwtabfile_tables):
    """Validate the values of all table sections.

    :param mwtabfile: Instance of :class:`~mwtab.mwtab.MWTabFile`.
    :type mwtabfile: :class:`~mwtab.mwtab.MWTabFile` or
                     :py:class:`collections.OrderedDict`
    :param data_section_key: Section key (either MS_METABOLITE_DATA, NMR_METABOLITE_DATA, or NMR_BINNED_DATA)
    :type data_section_key: :py:class:`str`
    :param mwtabfile_tables: Dictionary where the keys are table names and the values are the tables as pandas DataFrames.
    :type mwtabfile_tables: :py:class:`dict`
    """
    message_strings = {
        'Data': "[\"" + data_section_key + "\"][\"Data\"] \ METABOLITE_DATA" if not 'BINNED' in data_section_key else "[\"" + data_section_key + "\"][\"Data\"] \ BINNED_DATA",
        'Metabolites': "[\"" + data_section_key + "\"][\"Metabolites\"] \ METABOLITES",
        'Extended': "[\"" + data_section_key + "\"][\"Extended\"] \ EXTENDED_METABOLITE_DATA"
        }
    errors = []
    for table_name in mwtabfile.table_names:
        if table_name in mwtabfile[data_section_key]:
            records = mwtabfile[data_section_key][table_name]
            if len(records) > 0:
                headers = list(records[0].keys())
                if not all([list(data_dict.keys()) == headers for data_dict in records]):
                    # TODO test
                    message = ('Error: The table in the \"{}\" section does not have the same columns for every row. '
                               'The table, represented as a list of dictionaries, has dictionaries '
                               'with different keys. All dictionaries in the list should have the same keys.'.format(message_strings[table_name]))
                    errors.append(message)
            
            data_df = mwtabfile_tables[table_name]
            # temp_list = [duplicates_dict._JSON_DUPLICATE_KEYS__Jobj for duplicates_dict in mwtabfile[data_section_key][table_name]]
            # data_df = pandas.DataFrame.from_records(temp_list)
            # data_df = pandas.DataFrame.from_records(mwtabfile[data_section_key][table_name])
            
            # Look for completely null columns.
            null_columns = data_df.isna().all() | (data_df == "").all()
            null_columns = null_columns[null_columns]
            if len(null_columns) > 0:
                for column in null_columns.index:
                    message = "Warning: The column, \"{}\", in the {} section has all null values.".format(column, message_strings[table_name])
                    errors.append(message)
            
            # Look for overbalanced values, so if 90% of a column is dominated by a single value print a warning.
            for i, column in enumerate([column for column in data_df.columns if column != 'Metabolite']):
                temp_column = data_df.loc[:, column].astype(str)
                value_counts = temp_column.value_counts(dropna=False)
                value_counts = value_counts / value_counts.sum()
                if any(value_counts > .9) and len(value_counts) > 1:
                    message = ("Warning: The column, \"{}\", in the {} section may have incorrect values. "
                              "90% or more of the values are the same, but 10% or less are different.".format(column, message_strings[table_name]))
                    errors.append(message)
            
            # Look for duplicate rows.
            #TODO test.
            if data_df.duplicated().any():
                message = "Warning: There are duplicate rows in the {} section.".format(message_strings[table_name])
                errors.append(message)
    return errors


# TODO check the values of ion_mode and talk to hunter about restricting this to only pos or neg.
# Add docstring.
# Not sure there won't be MS that could have other modes.
def validate_polarity(mwtabfile, mwtabfile_tables):
    """
    
    :param mwtabfile_tables: Dictionary where the keys are table names and the values are the tables as pandas DataFrames.
    :type mwtabfile_tables: :py:class:`dict`
    """
    common_message = ('A single mwTab file is supposed to be restricted to a single analysis. '
                      'This means multiple MS runs under different settings should each be '
                      'in their own file.')
    errors = []
    if 'MS' in mwtabfile and 'ION_MODE' in mwtabfile['MS']:
        ion_mode = mwtabfile['MS']['ION_MODE']
        if ion_mode not in ['pos', 'neg', 'positive', 'negative']:
            errors.append('Error: The indicated ION_MODE should be either POSITIVE or NEGATIVE. ' + common_message)
    
    df = mwtabfile_tables['Metabolites']
    column_finder = column_finders['polarity']
    columns = {column:column.lower().strip() for column in df.columns}
    if column_matches := column_finder.name_dict_match(columns):
        for column_match in column_matches:
            pos_values = df.loc[:, column_match].apply(lambda x: 'pos' in x or '+' in x)
            neg_values = df.loc[:, column_match].apply(lambda x: 'neg' in x or '-' in x)
            if pos_values.any() and neg_values.any():
                errors.append(f'Error: The "{column_match}" column in the METABOLITES table '
                              'indicates multiple polarities in a single analysis, and '
                              'this should not be. ' + common_message)
                break
    return errors


def validate_file(mwtabfile, section_schema_mapping=section_schema_mapping, verbose=False, metabolites=True):
    """Validate ``mwTab`` formatted file.

    :param mwtabfile: Instance of :class:`~mwtab.mwtab.MWTabFile`.
    :type mwtabfile: :class:`~mwtab.mwtab.MWTabFile` or
                     :py:class:`collections.OrderedDict`
    :param dict section_schema_mapping: Dictionary that provides mapping between section name and schema definition.
    :param bool verbose: whether to be verbose or not.
    :param bool metabolites: whether to validate metabolites section.
    :return: Validated file and errors if verbose is False.
    :rtype: :py:class:`collections.OrderedDict`, _io.StringIO
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
    
    
    
    
    
    
    
    dict_for_Schema = OrderedDict()
    for section_key, section in mwtabfile.items():
        dict_for_Schema[section_key] = section
    try:
        base_schema.validate(dict_for_Schema)
    except Exception as e:
        messages = []
        for auto in e.autos:
            if auto.startswith("Wrong keys"):
                bad_key = match(r"Wrong keys ('.*') in .*", auto).group(1)
                auto = "Unknown or malformed keys " + bad_key + "."
            elif auto.startswith("Wrong key"):
                bad_key = match(r"Wrong key ('.*') .*", auto).group(1)
                auto = "Unknown or malformed key " + bad_key + "."
            messages.append(auto)
        errors.append("SCHEMA: There is an issue with the top level of the data.\n" + 
                      "\n".join(messages))
        # errors.append("SCHEMA: There is an issue with the top level of the data.\n" + 
        #               "\n".join([message for message in e.autos[1:] if message]))
    # Need to manually check this because Schema can't do conditional checks.
    if 'MS' in mwtabfile and 'CHROMATOGRAPHY' not in mwtabfile:
        if errors[-1].startswith("SCHEMA: There is an issue with the top level of the data."):
            errors[-1] = errors[-1] + "\nMissing key: 'CHROMATOGRAPHY'"
        else:
            errors.append("SCHEMA: There is an issue with the top level of the data.\nMissing key: 'CHROMATOGRAPHY'")

    # validate PROJECT, STUDY, ANALYSIS... and Schemas
    # TODO see if this catches multiple errors in a single section. For example, bad SUBJECT_TYPE and bad SUBJECT_SPECIES in the SUBJECT section.
    mwtabfile = {'MS': {'Data':2, 'Units':3}}
    errors = []
    for section_key, section in mwtabfile.items():
        try:
            schema = section_schema_mapping[section_key]
            # section = validate_section_schema(section, schema, section_key, error_stout)
            section, schema_errors = validate_section_schema(section, schema, section_key)
            errors.extend(schema_errors)
            # mwtabfile[section_key] = section
        except Exception as e:
            errors.append("SCHEMA: Section \"{}\" does not match the allowed schema. ".format(section_key) + str(e))






    # validate SUBJECT_SAMPLE_FACTORS
    # validate_subject_samples_factors(mwtabfile, error_stout)
    errors.extend(validate_subject_samples_factors(mwtabfile))
    errors.extend(validate_factors(mwtabfile))

    # validate ..._DATA sections
    data_section_key = mwtabfile.data_section_key
    if data_section_key:
        data_section_key = data_section_key[0]
        # validate_data(mwtabfile, data_section_key, error_stout, False)
        errors.extend(validate_data(mwtabfile, data_section_key, False, metabolites))

        if data_section_key in ("MS_METABOLITE_DATA", "NMR_METABOLITE_DATA"):
            # temp for testing
            if metabolites:
                if "Metabolites" in mwtabfile[data_section_key].keys():
                    errors.extend(validate_metabolites(mwtabfile, data_section_key))
                else:
                    errors.append("DATA: Missing METABOLITES section.")
        if "Extended" in mwtabfile[data_section_key].keys():
            errors.extend(validate_extended(mwtabfile, data_section_key))
        
        errors.extend(validate_samples(mwtabfile, data_section_key))
        errors.extend(validate_metabolite_names(mwtabfile, data_section_key))
        errors.extend(validate_table_values(mwtabfile, data_section_key))

    else:
        if "MS" in mwtabfile.keys():
            if not mwtabfile["MS"].get("MS_RESULTS_FILE"):
                errors.append("DATA: Missing MS_METABOLITE_DATA section or MS_RESULTS_FILE item in MS section.")
        elif "NM" in mwtabfile.keys():
            if not mwtabfile['NM'].get('NMR_RESULTS_FILE'):
                errors.append("DATA: Missing either NMR_METABOLITE_DATA or NMR_BINNED_DATA section or NMR_RESULTS_FILE item in NM section.")
    
    errors.extend(validate_metabolite_headers(mwtabfile))
    errors.extend(validate_header_lengths(mwtabfile))
    errors.extend(validate_sub_section_uniqueness(mwtabfile))

    # finish writing validation/error log
    if errors:
        print("Status: Contains Validation Errors", file=error_stout)
        print("Number Errors: {}\n".format(len(errors)), file=error_stout)
        try:
            print("Error Log:\n" + "\n".join(errors), file=error_stout)
        except UnicodeEncodeError:
            try:
                print("An error occurred when trying to print the error log with unicode. Trying UTF-8.")
                print(("Error Log:\n" + "\n".join(errors)).encode('utf-8'), file=error_stout)
            except Exception as e:
                print("An error occurred when trying to print the error log, so it could not be printed.")
                print("Error Log Exception:\n")
                traceback.print_exception(e, file=error_stout)
    else:
        print("Status: Passing", file=error_stout)

    if verbose:
        return mwtabfile, None
    else:
        return mwtabfile, error_stout.getvalue()

# TODO add checks for METABOLITES columns and try to give warnings for bad values, for example 'kegg_id' column should all be C00000, formula, inchi key,
# Change some error messages to only print once if the file is a tab file and not a JSON file. For example if a column name is off.

# When certain columns are found in METABOLITES, look for the implied pair and warn if it isn't there. For example, other_id and other_id_type  and retention_index and retention_index_type.
# Also look at the values in each pair and make sure they match, for instance AN000645 has values in other_id_type, but none in other_id.
# I put these in validate_metabolites.

# Look for ID values such as HMDB in the other_id column and suggest to make specific columns for them instead of putting them in other_id. Done
# Warn if a column name matches 2 different regexes in METABOLITES. 'retention time_m/z' in AN002889   AN004492 Feature@RT  Done

# Validate some values? Talk to Hunter. AN003426 has MS_METABOLITE_DATA:UNITS value as "counts" which seems wrong. Clearly peak area.

# Detect if multiple polarities in the same dataset and warn about it. UNSPECIFIED  MS:ION_MODE  or polarity column with both pos and neg

# Think about extending METABOLITES and EXTENDED blocks with an "Attributes" line like "Factors" in DATA block as a way to add more information about the columns themselves.

