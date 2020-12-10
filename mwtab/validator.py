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


def validate_subject_samples_factors(mwtabfile):
    """Validate ``SUBJECT_SAMPLE_FACTORS`` section.

    :param mwtabfile: Instance of :class:`~mwtab.mwtab.MWTabFile`.
    :type mwtabfile: :class:`~mwtab.mwtab.MWTabFile` or
                     :py:class:`collections.OrderedDict`
    """
    # TODO: REMOVE
    errors = list()

    for index, subject_sample_factor in enumerate(mwtabfile["SUBJECT_SAMPLE_FACTORS"]):
        if not subject_sample_factor["Subject ID"]:
            errors.append("SUBJECT_SAMPLE_FACTORS: Entry #{} missing Subject ID.".format(index+1))
        if not subject_sample_factor["Sample ID"]:
            errors.append("SUBJECT_SAMPLE_FACTORS: Entry #{} missing Sample ID.".format(index+1))
        if subject_sample_factor["Factors"]:
            for factor_key in subject_sample_factor["Factors"]:
                if not subject_sample_factor["Factors"][factor_key]:
                    errors.append("SUBJECT_SAMPLE_FACTORS: Entry #{} missing value for Factor {}.".format(index + 1, factor_key))
        if subject_sample_factor["Additional sample data"]:
            for additional_key in subject_sample_factor["Additional sample data"]:
                if not subject_sample_factor["Additional sample data"][additional_key]:
                    errors.append("SUBJECT_SAMPLE_FACTORS: Entry #{} missing value for Additional sample data {}.".format(index + 1, additional_key))

    return errors


def validate_data(mwtabfile, data_section_key):
    """Validates ``MS_METABOLITE_DATA``, ``NMR_METABOLITE_DATA``, and ``NMR_BINNED_DATA`` sections.

    :param mwtabfile: Instance of :class:`~mwtab.mwtab.MWTabFile`.
    :type mwtabfile: :class:`~mwtab.mwtab.MWTabFile` or
                     :py:class:`collections.OrderedDict`
    :param data_section_key: Section key (either MS_METABOLITE_DATA, NMR_METABOLITE_DATA, or NMR_BINNED_DATA)
    :type data_section_key: :py:class:`str`
    """
    # TODO: REMOVE
    errors = list()

    sample_id_set = {subject_sample_factor["Sample ID"] for subject_sample_factor in mwtabfile["SUBJECT_SAMPLE_FACTORS"]}

    for index, metabolite in enumerate(mwtabfile[data_section_key]["Data"]):
        if set(list(metabolite.keys())[1:]) != sample_id_set:
            errors.append("{} 1: Data entry #{} contains Sample IDs not matching those in SUBJECT_SAMPLE_FACTORS section.".format(data_section_key, index + 1))
        for data_point_key in metabolite.keys():
            if data_point_key != "Metabolite":
                try:
                    float(metabolite[data_point_key])
                except ValueError as e:
                    metabolite[data_point_key] = ""
                    errors.append("{} 2: Data entry #{} contains non-numeric value converted to \"\".".format(data_section_key, index + 1))

    return errors


def validate_metabolites(mwtabfile, data_section_key):
    """Validate ``METABOLITES`` section.

    :param mwtabfile: Instance of :class:`~mwtab.mwtab.MWTabFile`.
    :type mwtabfile: :class:`~mwtab.mwtab.MWTabFile` or
                     :py:class:`collections.OrderedDict`
    :param data_section_key: Section key (either MS_METABOLITE_DATA, NMR_METABOLITE_DATA, or NMR_BINNED_DATA)
    :type data_section_key: :py:class:`str`
    """
    # TODO: REMOVE
    errors = list()

    for index, metabolite in enumerate(mwtabfile[data_section_key]["Metabolites"]):
        for field_key in list(metabolite.keys())[1:]:
            if not any(k == field_key for k in METABOLITES_REGEXS.keys()):
                for regex_key in METABOLITES_REGEXS.keys():
                    if any(match(p, field_key) for p in METABOLITES_REGEXS[regex_key]):
                        errors.append("METABOLITES: Data entry #{} contains field name \"{}\" which matches a commonly used field name \"{}\".".format(index + 1, field_key, regex_key))
                        field_key = regex_key
                        break

    return errors


def validate_extended(mwtabfile, data_section_key):
    """Validate ``EXTENDED_MS_METABOLITE_DATA``, ``EXTENDED_NMR_METABOLITE_DATA``, and ``EXTENDED_NMR_BINNED_DATA`` sections.

    :param mwtabfile: Instance of :class:`~mwtab.mwtab.MWTabFile`.
    :type mwtabfile: :class:`~mwtab.mwtab.MWTabFile` or
                     :py:class:`collections.OrderedDict`
    :param data_section_key: Section key (either MS_METABOLITE_DATA, NMR_METABOLITE_DATA, or NMR_BINNED_DATA)
    :type data_section_key: :py:class:`str`
    """
    # TODO: REMOVE
    errors = list()

    sample_id_set = {subject_sample_factor["Sample ID"] for subject_sample_factor in
                     mwtabfile["SUBJECT_SAMPLE_FACTORS"]}

    for index, extended_data in enumerate(mwtabfile[data_section_key]["Extended"]):
        if "sample_id" not in extended_data.keys():
            errors.append("EXTENDED_{}: Data entry #{} missing Sample ID.".format(data_section_key, index + 1))
        elif not extended_data["sample_id"] in sample_id_set:
            errors.append("EXTENDED_{}: Data entry #{} contains Sample ID \"{}\" not found in SUBJECT_SAMPLE_FACTORS section.".format(data_section_key, index + 1, extended_data["sample_id"]))

    return errors


def validate_section_schema(section, schema, section_key):
    """Validate section of ``mwTab`` formatted file.

    :param section: Section of :class:`~mwtab.mwtab.MWTabFile`.
    :type section: :py:class:`collections.OrderedDict`
    :param schema: Schema definition.
    :type schema:
    :param section_key: Section key.
    :type section_key: :py:class:`str`

    :return: Validated section.
    :rtype: :py:class:`collections.OrderedDict`
    """
    # TODO: REMOVE
    errors = list()

    if section_key in ITEM_SECTIONS:
        for key in section.keys():
            if not section[key]:
                errors.append("{}: Contains item \"{}\" with null value.".format(section_key, key))
                del section[key]

    return schema.validate(section), errors


def validate_file(mwtabfile, section_schema_mapping=section_schema_mapping):
    """Validate ``mwTab`` formatted file.

    :param mwtabfile: Instance of :class:`~mwtab.mwtab.MWTabFile`.
    :type mwtabfile: :class:`~mwtab.mwtab.MWTabFile` or
                     :py:class:`collections.OrderedDict`
    :param section_schema_mapping: Dictionary that provides mapping between section name and schema definition.
`    :type section_schema_mapping: :py:class:`dict`
`
    :return: Validated file.
    :rtype: :py:class:`collections.OrderedDict`
    """
    errors = list()
    validated_mwtabfile = deepcopy(OrderedDict(mwtabfile))

    for section_key, section in mwtabfile.items():
        try:
            schema = section_schema_mapping[section_key]
            validate_section_schema(section, schema, section_key)
            validated_mwtabfile[section_key] = section
        except Exception as e:
            errors.append("SCHEMA: Section \"{}\" does not match the allowed schema. ".format(section_key) + str(e))

    errors.extend(validate_subject_samples_factors(validated_mwtabfile))

    data_section_key = list(set(validated_mwtabfile.keys()) & {"MS_METABOLITE_DATA", "NMR_METABOLITE_DATA", "NMR_BINNED_DATA"})
    if data_section_key:
        data_section_key = data_section_key[0]
        errors.extend(validate_data(validated_mwtabfile, data_section_key))
        if data_section_key in ("MS_METABOLITE_DATA", "NMR_METABOLITE_DATA"):
            if "Metabolites" in validated_mwtabfile[data_section_key].keys():
                errors.extend(validate_metabolites(validated_mwtabfile, data_section_key))
            else:
                errors.append("DATA: Missing METABOLITES section.")
        if "Extended" in validated_mwtabfile[data_section_key].keys():
            errors.extend(validate_extended(validated_mwtabfile, data_section_key))

    else:
        if "MS" in validated_mwtabfile.keys():
            errors.append("DATA: Missing MS_METABOLITE_DATA section.")
        elif "NM":
            errors.append("DATA: Missing either NMR_METABOLITE_DATA or NMR_BINNED_DATA section.")

    return validated_mwtabfile, errors
