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
from .mwschema import section_schema_mapping
from re import match
import io
import sys
import logging
import warnings


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


def validate_subject_samples_factors(mwtabfile, error_stream):
    """Validate ``SUBJECT_SAMPLE_FACTORS`` section.

    :param mwtabfile: Instance of :class:`~mwtab.mwtab.MWTabFile`.
    :type mwtabfile: :class:`~mwtab.mwtab.MWTabFile` or
                     :py:class:`collections.OrderedDict`
    :param error_stream: IO stream to print errors to.
    """
    for index, subject_sample_factor in enumerate(mwtabfile["SUBJECT_SAMPLE_FACTORS"]):
        if not subject_sample_factor["Subject ID"]:
            print("SUBJECT_SAMPLE_FACTORS: Entry #{} missing Subject ID.".format(index+1), file=error_stream)
        if not subject_sample_factor["Sample ID"]:
            print("SUBJECT_SAMPLE_FACTORS: Entry #{} missing Sample ID.".format(index + 1), file=error_stream)
        if subject_sample_factor.get("Factors"):
            for factor_key in subject_sample_factor["Factors"]:
                if not subject_sample_factor["Factors"][factor_key]:
                    print(
                        "SUBJECT_SAMPLE_FACTORS: Entry #{} missing value for Factor {}.".format(index + 1, factor_key),
                        file=error_stream
                    )
        if subject_sample_factor.get("Additional sample data"):
            for additional_key in subject_sample_factor["Additional sample data"]:
                if not subject_sample_factor["Additional sample data"][additional_key]:
                    print(
                        "SUBJECT_SAMPLE_FACTORS: Entry #{} missing value for Additional sample data {}.".format(
                            index + 1, additional_key
                        ),
                        file=error_stream
                    )


def validate_data(mwtabfile, data_section_key, error_stream, null_values):
    """Validates ``MS_METABOLITE_DATA``, ``NMR_METABOLITE_DATA``, and ``NMR_BINNED_DATA`` sections.

    :param mwtabfile: Instance of :class:`~mwtab.mwtab.MWTabFile`.
    :type mwtabfile: :class:`~mwtab.mwtab.MWTabFile` or
                     :py:class:`collections.OrderedDict`
    :param data_section_key: Section key (either MS_METABOLITE_DATA, NMR_METABOLITE_DATA, or NMR_BINNED_DATA)
    :type data_section_key: :py:class:`str`
    :param error_stream: IO stream to print errors to.
    :param bool null_values: whether null values are present.
    """
    sample_id_set = {subject_sample_factor["Sample ID"] for subject_sample_factor in mwtabfile["SUBJECT_SAMPLE_FACTORS"]}

    for metabolite in mwtabfile[data_section_key]["Data"]:
        if set(list(metabolite.keys())[1:]) != sample_id_set:
            print("{}: Metabolite/Bin Range \"{}\" missing data entry for {} samples".format(data_section_key, metabolite[list(metabolite.keys())[0]], len(sample_id_set.symmetric_difference(set(list(metabolite.keys())[1:])))), file=error_stream)
        if null_values:
            for data_point_key in metabolite.keys():
                if data_point_key != "Metabolite":
                    try:
                        float(metabolite[data_point_key])
                    except ValueError as e:
                        metabolite[data_point_key] = ""
                        print("{}: Metabolite/Bin Range \"{}\" contains non-numeric value converted to \"\".".format(data_section_key, metabolite[list(metabolite.keys())[0]]), file=error_stream)


def validate_metabolites(mwtabfile, data_section_key, error_stream):
    """Validate ``METABOLITES`` section.

    :param mwtabfile: Instance of :class:`~mwtab.mwtab.MWTabFile`.
    :type mwtabfile: :class:`~mwtab.mwtab.MWTabFile` or
                     :py:class:`collections.OrderedDict`
    :param data_section_key: Section key (either MS_METABOLITE_DATA, NMR_METABOLITE_DATA, or NMR_BINNED_DATA)
    :type data_section_key: :py:class:`str`
    :param error_stream: IO stream to print errors to.
    """
    for index, metabolite in enumerate(mwtabfile[data_section_key]["Metabolites"]):
        for field_key in list(metabolite.keys())[1:]:
            if not any(k == field_key for k in METABOLITES_REGEXS.keys()):
                for regex_key in METABOLITES_REGEXS.keys():
                    if any(match(p, field_key) for p in METABOLITES_REGEXS[regex_key]):
                        print("METABOLITES: Data entry #{} contains field name \"{}\" which matches a commonly used field name \"{}\".".format(index + 1, field_key, regex_key), file=error_stream)
                        field_key = regex_key
                        break


def validate_extended(mwtabfile, data_section_key, error_stream):
    """Validate ``EXTENDED_MS_METABOLITE_DATA``, ``EXTENDED_NMR_METABOLITE_DATA``, and ``EXTENDED_NMR_BINNED_DATA`` sections.

    :param mwtabfile: Instance of :class:`~mwtab.mwtab.MWTabFile`.
    :type mwtabfile: :class:`~mwtab.mwtab.MWTabFile` or
                     :py:class:`collections.OrderedDict`
    :param data_section_key: Section key (either MS_METABOLITE_DATA, NMR_METABOLITE_DATA, or NMR_BINNED_DATA)
    :type data_section_key: :py:class:`str`
    :param error_stream: IO stream to print errors to.
    """
    sample_id_set = {subject_sample_factor["Sample ID"] for subject_sample_factor in
                     mwtabfile["SUBJECT_SAMPLE_FACTORS"]}

    for index, extended_data in enumerate(mwtabfile[data_section_key]["Extended"]):
        if "sample_id" not in extended_data.keys():
            print("EXTENDED_{}: Data entry #{} missing Sample ID.".format(data_section_key, index + 1),
                  file=error_stream)
        elif not extended_data["sample_id"] in sample_id_set:
            print(
                "EXTENDED_{}: Data entry #{} contains Sample ID \"{}\" not found in SUBJECT_SAMPLE_FACTORS section.".format(
                    data_section_key, index + 1, extended_data["sample_id"]
                ),
                file=error_stream
            )


def validate_section_schema(section, schema, section_key, error_stream):
    """Validate section of ``mwTab`` formatted file.

    :param section: Section of :class:`~mwtab.mwtab.MWTabFile`.
    :type section: :py:class:`collections.OrderedDict`
    :param schema: Schema definition.
    :type schema: :py:class:`~schema.schema`
    :param str section_key: Section key.
    :param error_stream: IO stream to print errors to.
    :return: Validated section.
    :rtype: :py:class:`collections.OrderedDict`
    """
    if section_key in ITEM_SECTIONS:
        for key in section.keys():
            if not section[key]:
                print("{}: Contains item \"{}\" with null value.".format(section_key, key), file=error_stream)
                del section[key]

    return schema.validate(section)


def validate_file(mwtabfile, section_schema_mapping=section_schema_mapping, verbose=False, metabolites=True):
    """Validate ``mwTab`` formatted file.

    :param mwtabfile: Instance of :class:`~mwtab.mwtab.MWTabFile`.
    :type mwtabfile: :class:`~mwtab.mwtab.MWTabFile` or
                     :py:class:`collections.OrderedDict`
    :param dict section_schema_mapping: Dictionary that provides mapping between section name and schema definition.
    :param bool verbose: whether to be verbose or not.
    :param bool metabolites: whether to validate metabolites section.
    :return: Validated file.
    :rtype: :py:class:`collections.OrderedDict`
    """
    # setup
    if not verbose:
        error_stout = io.StringIO()
    else:
        error_stout = sys.stdout
    validated_mwtabfile = deepcopy(OrderedDict(mwtabfile))

    # validate PROJECT, STUDY, ANALYSIS... and Schemas
    for section_key, section in mwtabfile.items():
        try:
            schema = section_schema_mapping[section_key]
            section = validate_section_schema(section, schema, section_key, error_stout)
            validated_mwtabfile[section_key] = section
        except Exception as e:
            print("SCHEMA: Section \"{}\" does not match the allowed schema. ".format(section_key) + str(e), file=error_stout)

    # validate SUBJECT_SAMPLE_FACTORS
    validate_subject_samples_factors(validated_mwtabfile, error_stout)

    data_section_key = list(set(validated_mwtabfile.keys()) & {"MS_METABOLITE_DATA", "NMR_METABOLITE_DATA", "NMR_BINNED_DATA"})
    if data_section_key:
        data_section_key = data_section_key[0]
        validate_data(validated_mwtabfile, data_section_key, error_stout, False)
        if data_section_key in ("MS_METABOLITE_DATA", "NMR_METABOLITE_DATA"):
            # temp for testing
            if metabolites:
                if "Metabolites" in validated_mwtabfile[data_section_key].keys():
                    validate_metabolites(validated_mwtabfile, data_section_key, error_stout)
                else:
                    print("DATA: Missing METABOLITES section.", file=error_stout)
        if "Extended" in validated_mwtabfile[data_section_key].keys():
            validate_extended(validated_mwtabfile, data_section_key, error_stout)

    else:
        if "MS" in validated_mwtabfile.keys():
            if not validated_mwtabfile["MS"].get("MS_RESULTS_FILE"):
                print("DATA: Missing MS_METABOLITE_DATA section or MS_RESULTS_FILE item in MS section.", file=error_stout)
        elif "NM":
            print("DATA: Missing either NMR_METABOLITE_DATA or NMR_BINNED_DATA section.", file=error_stout)

    errors = ""
    if not verbose:
        errors = error_stout.getvalue()

    return validated_mwtabfile, errors
