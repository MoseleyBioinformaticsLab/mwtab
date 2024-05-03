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

from copy import deepcopy
from collections import OrderedDict
from datetime import datetime
from re import match
import io
import sys
import traceback

from .mwschema import section_schema_mapping, base_schema, _duplicate_key_list

import mwtab


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


def validate_subject_samples_factors(mwtabfile):
    """Validate ``SUBJECT_SAMPLE_FACTORS`` section.

    :param mwtabfile: Instance of :class:`~mwtab.mwtab.MWTabFile`.
    :type mwtabfile: :class:`~mwtab.mwtab.MWTabFile` or
                     :py:class:`collections.OrderedDict`
    """
    subject_samples_factors_errors = list()

    for index, subject_sample_factor in enumerate(mwtabfile["SUBJECT_SAMPLE_FACTORS"]):
        if not subject_sample_factor["Subject ID"]:
            subject_samples_factors_errors.append(
                "SUBJECT_SAMPLE_FACTORS: Entry #{} missing Subject ID.".format(index+1)
            )
        if not subject_sample_factor["Sample ID"]:
            subject_samples_factors_errors.append(
                "SUBJECT_SAMPLE_FACTORS: Entry #{} missing Sample ID.".format(index + 1)
            )
        if subject_sample_factor.get("Factors"):
            for factor_key in subject_sample_factor["Factors"]:
                if not subject_sample_factor["Factors"][factor_key]:
                    subject_samples_factors_errors.append(
                        "SUBJECT_SAMPLE_FACTORS: Entry #{} missing value for Factor {}.".format(index + 1, factor_key)
                    )
        if subject_sample_factor.get("Additional sample data"):
            for additional_key in subject_sample_factor["Additional sample data"]:
                if not subject_sample_factor["Additional sample data"][additional_key]:
                    subject_samples_factors_errors.append(
                        "SUBJECT_SAMPLE_FACTORS: Entry #{} missing value for Additional sample data {}.".format(
                            index + 1, additional_key
                        )
                    )
            
            duplicate_keys = [key for key, value in subject_sample_factor["Additional sample data"].items() 
                              if isinstance(value, _duplicate_key_list)]
            if duplicate_keys:
                subject_samples_factors_errors.append("SUBJECT_SAMPLE_FACTORS: Entry #" + str(index + 1) + 
                                                      " has the following duplicate keys:\n\t" + 
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
    if mwtabfile[data_section_key]["Data"]:
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
        data_errors.append("SUBJECT_SAMPLE_FACTORS: Section missing sample ID(s) {} found in {} section.".format(
            data_sample_id_set - subject_sample_factors_sample_id_set,
            data_section_key
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
        
        # Check whther the metabolite is in the Metabolites section.
        if metabolites and "Metabolites" in mwtabfile[data_section_key]:
            metabolite = metabolite_dict["Metabolite"]
            if metabolite not in metabolites_in_metabolites_section:
                data_errors.append("DATA: Data entry #{}, \"{}\", is not in the Metabolites section.".format(index + 1, metabolite))
        
        # Check if there are null values.
        if null_values:
            for data_point_key in metabolite_dict.keys():
                if data_point_key != "Metabolite":
                    try:
                        float(metabolite_dict[data_point_key])
                    except ValueError as e:
                        metabolite_dict[data_point_key] = ""
                        data_errors.append(
                            "{}: Data entry #{} contains non-numeric value converted to \"\".".format(data_section_key, index + 1))

    return data_errors


def validate_metabolites(mwtabfile, data_section_key):
    """Validate ``METABOLITES`` section.

    :param mwtabfile: Instance of :class:`~mwtab.mwtab.MWTabFile`.
    :type mwtabfile: :class:`~mwtab.mwtab.MWTabFile` or
                     :py:class:`collections.OrderedDict`
    :param data_section_key: Section key (either MS_METABOLITE_DATA, NMR_METABOLITE_DATA, or NMR_BINNED_DATA)
    :type data_section_key: :py:class:`str`
    """
    metabolites_errors = list()
    
    metabolites_in_data_section = [data_dict["Metabolite"] for data_dict in mwtabfile[data_section_key]["Data"]]

    for index, metabolite_dict in enumerate(mwtabfile[data_section_key]["Metabolites"]):
        # Check whether the metabolite is in the Data section.
        metabolite = metabolite_dict["Metabolite"]
        if metabolite not in metabolites_in_data_section:
            metabolites_errors.append("METABOLITES: Data entry #{}, \"{}\", is not in the Data section.".format(index + 1, metabolite))
        
        # Check if fields are recognized variations and report the standardized name to the user.
        for field_key in list(metabolite_dict.keys())[1:]:
            if not any(k == field_key for k in METABOLITES_REGEXS.keys()):
                for regex_key in METABOLITES_REGEXS.keys():
                    if any(match(p, field_key) for p in METABOLITES_REGEXS[regex_key]):
                        metabolites_errors.append("METABOLITES: Data entry #{} contains field name \"{}\" which matches a commonly used field name \"{}\".".format(index + 1, field_key, regex_key))
                        field_key = regex_key
                        break

    return metabolites_errors


def validate_extended(mwtabfile, data_section_key):
    """Validate ``EXTENDED_MS_METABOLITE_DATA``, ``EXTENDED_NMR_METABOLITE_DATA``, and ``EXTENDED_NMR_BINNED_DATA`` sections.

    :param mwtabfile: Instance of :class:`~mwtab.mwtab.MWTabFile`.
    :type mwtabfile: :class:`~mwtab.mwtab.MWTabFile` or
                     :py:class:`collections.OrderedDict`
    :param data_section_key: Section key (either MS_METABOLITE_DATA, NMR_METABOLITE_DATA, or NMR_BINNED_DATA)
    :type data_section_key: :py:class:`str`
    """
    extended_errors = list()

    sample_id_set = {subject_sample_factor["Sample ID"] for subject_sample_factor in
                     mwtabfile["SUBJECT_SAMPLE_FACTORS"]}

    for index, extended_data in enumerate(mwtabfile[data_section_key]["Extended"]):
        if "sample_id" not in extended_data.keys():
            extended_errors.append("EXTENDED_{}: Data entry #{} missing Sample ID.".format(data_section_key, index + 1))
        elif not extended_data["sample_id"] in sample_id_set:
            extended_errors.append(
                "EXTENDED_{}: Data entry #{} contains Sample ID \"{}\" not found in SUBJECT_SAMPLE_FACTORS section.".format(
                    data_section_key, index + 1, extended_data["sample_id"]
                ))

    return extended_errors


def validate_section_schema(section, schema, section_key, cleaning=False):
    """Validate section of ``mwTab`` formatted file.
    
    :param section: Section of :class:`~mwtab.mwtab.MWTabFile`.
    :type section: :py:class:`collections.OrderedDict`
    :param schema: Schema definition.
    :type schema: :py:class:`~schema.schema`
    :param str section_key: Section key.
    :return: Validated section.
    :rtype: :py:class:`collections.OrderedDict`
    """
    schema_errors = list()

    keys_to_delete = []
    if section_key in ITEM_SECTIONS:
        for key in section.keys():
            if not section[key]:
                schema_errors.append("{}: Contains item \"{}\" with null value.".format(section_key, key))
                keys_to_delete.append(key)

    if cleaning:
        for key in keys_to_delete:
            del section[key]

    return schema.validate(section), schema_errors


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
    validated_mwtabfile = deepcopy(OrderedDict(mwtabfile))

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
    
    try:
        base_schema.validate(validated_mwtabfile)
    except Exception as e:
        errors.append("SCHEMA: There is an issue with the top level of the data." + 
                      "\n".join([message for message in e.autos[1:] if message]))

    # validate PROJECT, STUDY, ANALYSIS... and Schemas
    for section_key, section in mwtabfile.items():
        try:
            schema = section_schema_mapping[section_key]
            # section = validate_section_schema(section, schema, section_key, error_stout)
            section, schema_errors = validate_section_schema(section, schema, section_key)
            errors.extend(schema_errors)
            validated_mwtabfile[section_key] = section
        except Exception as e:
            errors.append("SCHEMA: Section \"{}\" does not match the allowed schema. ".format(section_key) + str(e))

    # validate SUBJECT_SAMPLE_FACTORS
    # validate_subject_samples_factors(validated_mwtabfile, error_stout)
    errors.extend(validate_subject_samples_factors(validated_mwtabfile))

    # validate ..._DATA sections
    data_section_key = list(set(validated_mwtabfile.keys()) &
                            {"MS_METABOLITE_DATA", "NMR_METABOLITE_DATA", "NMR_BINNED_DATA"})
    if data_section_key:
        data_section_key = data_section_key[0]
        # validate_data(validated_mwtabfile, data_section_key, error_stout, False)
        errors.extend(validate_data(validated_mwtabfile, data_section_key, False, metabolites))

        if data_section_key in ("MS_METABOLITE_DATA", "NMR_METABOLITE_DATA"):
            # temp for testing
            if metabolites:
                if "Metabolites" in validated_mwtabfile[data_section_key].keys():
                    errors.extend(validate_metabolites(validated_mwtabfile, data_section_key))
                else:
                    errors.append("DATA: Missing METABOLITES section.")
        if "Extended" in validated_mwtabfile[data_section_key].keys():
            errors.extend(validate_extended(validated_mwtabfile, data_section_key))

    else:
        if "MS" in validated_mwtabfile.keys():
            if not validated_mwtabfile["MS"].get("MS_RESULTS_FILE"):
                errors.append("DATA: Missing MS_METABOLITE_DATA section or MS_RESULTS_FILE item in MS section.")
        elif "NM" in validated_mwtabfile.keys():
            if not validated_mwtabfile['NM'].get('NMR_RESULTS_FILE'):
                errors.append("DATA: Missing either NMR_METABOLITE_DATA or NMR_BINNED_DATA section or NMR_RESULTS_FILE item in NM secction.")

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
        return validated_mwtabfile, None
    else:
        return validated_mwtabfile, error_stout.getvalue()
