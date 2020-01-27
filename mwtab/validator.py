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
from .mwschema import section_schema_mapping


def _validate_samples_factors(mwtabfile, validate_samples=True, validate_factors=True):
    """Validate ``Samples`` and ``Factors`` identifiers across the file.

    :param mwtabfile: Instance of :class:`~mwtab.mwtab.MWTabFile`.
    :type mwtabfile: :class:`~mwtab.mwtab.MWTabFile`
    :param validate_samples:
    :type validate_samples: :py:obj:`bool`
    :param validate_factors:
    :type validate_factors: :py:obj:`bool`
    :return: List of errors (["None"] if no errors).
    :rtype: :py:obj:`list`
    """
    errors = []

    from_subject_samples = {i["local_sample_id"] for i in mwtabfile["SUBJECT_SAMPLE_FACTORS"]["SUBJECT_SAMPLE_FACTORS"]}
    from_subject_factors = {i["factors"] for i in mwtabfile["SUBJECT_SAMPLE_FACTORS"]["SUBJECT_SAMPLE_FACTORS"]}

    if "" in from_subject_samples:
        errors.append(ValueError("Sample with no Sample ID (\"\") in `SUBJECT_SAMPLE_FACTOR` block."))
    if "" in from_subject_factors:
        errors.append(ValueError("Sample with no Factor(s) (\"\") in `SUBJECT_SAMPLE_FACTOR` block."))

    if validate_samples:

        if "MS_METABOLITE_DATA" in mwtabfile:
            try:
                from_metabolite_data_samples = set(mwtabfile["MS_METABOLITE_DATA"]["MS_METABOLITE_DATA_START"]["Samples"])
                if not from_metabolite_data_samples.issubset(from_subject_samples):
                    # TODO: Using list/OrderedDict to separate missing sample ids from tab errors
                    if "" in from_metabolite_data_samples:
                        errors.append(ValueError(
                            "Sample with no Sample ID (\"\") in `MS_METABOLITE_DATA` block (usually caused by "
                            "extraneous TAB at the end of line)."))
                    if any(val for val in from_metabolite_data_samples if val != ""):
                        error_msg = "`MS_METABOLITE_DATA` block contains additional samples not found in " \
                                    "`SUBJECT_SAMPLE_FACTORS` block.\n\tAdditional samples: {}".format(
                                        str(from_metabolite_data_samples.difference(from_subject_samples)))
                        errors.append(ValueError(error_msg))
            except KeyError:
                errors.append(KeyError("Missing key `Samples` in `MS_METABOLITE_DATA` block."))

        if "NMR_BINNED_DATA" in mwtabfile:
            try:
                from_nmr_binned_data_samples = set(mwtabfile["NMR_BINNED_DATA"]["NMR_BINNED_DATA_START"]["Fields"][1:])
                if not from_nmr_binned_data_samples.issubset(from_subject_samples):
                    if "" in from_nmr_binned_data_samples:
                        errors.append(ValueError(
                            "Sample with no sample ID (\"\") in `NMR_BINNED_DATA` block (usually caused by "
                            "extraneous TAB at the end of line)."))
                       if any(val for val in from_nmr_binned_data_samples if val != ""):
                        errors.append(ValueError(
                            "`NMR_BINNED_DATA` block contains additional samples not found in `SUBJECT_SAMPLE_FACTORS` "
                            "block.\n\tAdditional samples: {}".format(
                                from_nmr_binned_data_samples.difference(from_subject_samples))))
            except KeyError:
                errors.append(KeyError("Missing key `Bin range(ppm)` in `NMR_BINNED_DATA` block."))

    if validate_factors:

        if "MS_METABOLITE_DATA" in mwtabfile:
            try:
                from_metabolite_data_factors = set(mwtabfile["MS_METABOLITE_DATA"]["MS_METABOLITE_DATA_START"]["Factors"])
                if not from_metabolite_data_factors.issubset(from_subject_factors):
                    if "" in from_metabolite_data_factors:
                        errors.append(ValueError(
                            "Sample with no factors (\"\") in `MS_METABOLITE_DATA` block (usually caused by "
                            "extraneous TAB at the end of line)."))
                    if any(val for val in from_metabolite_data_factors if val != ""):
                        errors.append(ValueError(
                            "`MS_METABOLITE_DATA` block contains additional factors not found in "
                            "`SUBJECT_SAMPLE_FACTORS` block.\n\t Additional factors: {}".format(
                                from_metabolite_data_factors.difference(from_subject_samples))))
            except KeyError:
                errors.append(KeyError("Missing key `Factors` in `MS_METABOLITE_DATA` block."))

    return errors


def _validate_metabolites(mwtabfile):
    """Validate metabolite ``Features`` across the file.

    :param mwtabfile: Instance of :class:`~mwtab.mwtab.MWTabFile`.
    :type mwtabfile: :class:`~mwtab.mwtab.MWTabFile`
    :return: List of errors (["None"] if no errors).
    :rtype: :py:obj:`list`
    """
    errors = []

    try:
        from_metabolite_data_features = {i["metabolite_name"] for i in mwtabfile["MS_METABOLITE_DATA"]["MS_METABOLITE_DATA_START"]["DATA"]}
    except KeyError:
        errors.append(KeyError("Missing key `metabolite_name` in `MS_METABOLITE_DATA` block."))
    try:
        from_metabolites_features = {i["metabolite_name"] for i in mwtabfile["METABOLITES"]["METABOLITES_START"]["DATA"]}
    except KeyError:
        errors.append(KeyError("Missing key `metabolite_name` in `MS_METABOLITE_DATA` block."))

    if not errors:
        if "" in from_metabolite_data_features:
            errors.append(ValueError("Feature with no name (\"\") in `MS_METABOLITE_DATA` block."))
        if "" in from_metabolites_features:
            errors.append(ValueError("Feature with no name (\"\") in `METABOLITES` block."))

        if from_metabolite_data_features != from_metabolites_features:
            if any(val for val in from_metabolite_data_features if val != ""):
                errors.append(ValueError(
                    "`MS_METABOLITE_DATA` block contains additional features not found in `METABOLITES` block.\n\t"
                    "Additional features: {}".format(
                        from_metabolite_data_features.difference(from_metabolites_features))))

    return errors


def validate_section(section, schema):
    """Validate section of ``mwTab`` formatted file.

    :param section: Section of :class:`~mwtab.mwtab.MWTabFile`.
    :param schema: Schema definition.
    :return: Validated section.
    :rtype: :py:class:`collections.OrderedDict`
    """
    return schema.validate(section)


def validate_file(mwtabfile, section_schema_mapping=section_schema_mapping, validate_samples=True,
                  validate_factors=True, validate_features=True, validate_schema=True):
    """Validate entire ``mwTab`` formatted file one section at a time.

    :param mwtabfile: Instance of :class:`~mwtab.mwtab.MWTabFile`.
    :type mwtabfile: :class:`~mwtab.mwtab.MWTabFile`
    :param dict section_schema_mapping: Dictionary that provides mapping between section name and schema definition.
    :param validate_samples: Make sure that sample ids are consistent across file.
    :type validate_samples: :py:obj:`True` or :py:obj:`False`
    :param validate_factors: Make sure that factors are consistent across file.
    :type validate_factors: :py:obj:`True` or :py:obj:`False`
    :return: Validated file.
    :rtype: :py:class:`collections.OrderedDict`
    """
    errors = []

    if validate_samples or validate_factors:
        errors.extend(_validate_samples_factors(mwtabfile, validate_samples, validate_factors))

    if mwtabfile.get("METABOLITES") and validate_features:
        errors.extend(_validate_metabolites(mwtabfile))

    if validate_schema:
        for section_key, section in mwtabfile.items():
            try:
                schema = section_schema_mapping[section_key]
                validate_section(section=section, schema=schema)
            except Exception as e:
                errors.append(e)

    return errors
