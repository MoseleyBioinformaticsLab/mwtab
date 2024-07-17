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

import json_duplicate_keys as jdks

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
        factors_dict_1 = {i["Sample ID"]: i["Factors"] for i in mwtabfile["SUBJECT_SAMPLE_FACTORS"]}
        factors_dict_2 = {sample: factors for sample, factors in mwtabfile._factors.items()}
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
    """Validate that the samples in SUBJECT_SAMPLE_FACTORS and METABOLITE_DATA are the same.
    
    :param mwtabfile: Instance of :class:`~mwtab.mwtab.MWTabFile`.
    :type mwtabfile: :class:`~mwtab.mwtab.MWTabFile` or
                     :py:class:`collections.OrderedDict`
    :param data_section_key: Section key (either MS_METABOLITE_DATA, NMR_METABOLITE_DATA, or NMR_BINNED_DATA)
    :type data_section_key: :py:class:`str`
    """
    errors = []
    
    diff = None
    if "SUBJECT_SAMPLE_FACTORS" in mwtabfile:
        ssf_samples = {ssf_dict["Sample ID"] for ssf_dict in mwtabfile["SUBJECT_SAMPLE_FACTORS"] if 
                       "Sample ID" in ssf_dict and ssf_dict["Sample ID"]}
    else:
        return errors
    
    if mwtabfile._samples:
        diff = set(mwtabfile._samples) ^ ssf_samples
    elif "Data" in mwtabfile[data_section_key] and mwtabfile[data_section_key]["Data"]:
        samples = {sample for sample in mwtabfile[data_section_key]["Data"][0] if sample != "Metabolite"}
        diff = samples ^ ssf_samples
    if diff:
        errors.append(
            "Samples: The Sample ID's in the SUBJECT_SAMPLE_FACTORS do not match what is in the METABOLITE_DATA section."
        )
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
                    except ValueError as e:
                        # metabolite_dict[data_point_key] = ""
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
        errors.append(message.format(name, "[\"" + data_section_key + "\"][\"Data\"] \ METABOLITE_DATA"))
    for name in extended_section_bad_names:
        errors.append(message.format(name, "[\"" + data_section_key + "\"][\"Extended\"] \ EXTENDED_METABOLITE_DATA"))
    
    return errors
    


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
    validated_mwtabfile.source = mwtabfile.source
    validated_mwtabfile._factors = mwtabfile._factors
    validated_mwtabfile._samples = mwtabfile._samples
    validated_mwtabfile._metabolite_header = mwtabfile._metabolite_header
    validated_mwtabfile._extended_metabolite_header = mwtabfile._extended_metabolite_header
    validated_mwtabfile._binned_header = mwtabfile._binned_header
    validated_mwtabfile._short_headers = mwtabfile._short_headers
    validated_mwtabfile._duplicate_sub_sections = mwtabfile._duplicate_sub_sections

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
    errors.extend(validate_factors(validated_mwtabfile))

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
        
        errors.extend(validate_samples(validated_mwtabfile, data_section_key))
        
        errors.extend(validate_metabolite_names(validated_mwtabfile, data_section_key))

    else:
        if "MS" in validated_mwtabfile.keys():
            if not validated_mwtabfile["MS"].get("MS_RESULTS_FILE"):
                errors.append("DATA: Missing MS_METABOLITE_DATA section or MS_RESULTS_FILE item in MS section.")
        elif "NM" in validated_mwtabfile.keys():
            if not validated_mwtabfile['NM'].get('NMR_RESULTS_FILE'):
                errors.append("DATA: Missing either NMR_METABOLITE_DATA or NMR_BINNED_DATA section or NMR_RESULTS_FILE item in NM secction.")
    
    errors.extend(validate_metabolite_headers(validated_mwtabfile))
    errors.extend(validate_header_lengths(validated_mwtabfile))
    errors.extend(validate_sub_section_uniqueness(validated_mwtabfile))

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
