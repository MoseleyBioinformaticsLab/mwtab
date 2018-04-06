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
    :return: None
    :rtype: :py:obj:`None`
    """
    from_subject_samples = {i["local_sample_id"] for i in mwtabfile["SUBJECT_SAMPLE_FACTORS"]["SUBJECT_SAMPLE_FACTORS"]}
    from_subject_factors = {i["factors"] for i in mwtabfile["SUBJECT_SAMPLE_FACTORS"]["SUBJECT_SAMPLE_FACTORS"]}

    if validate_samples:

        if "MS_METABOLITE_DATA" in mwtabfile:
            from_metabolite_data_samples = set(mwtabfile["MS_METABOLITE_DATA"]["MS_METABOLITE_DATA_START"]["Samples"])
            assert from_subject_samples == from_metabolite_data_samples

        if "NMR_BINNED_DATA" in mwtabfile:
            from_nmr_binned_data_samples = set(mwtabfile["NMR_BINNED_DATA"]["NMR_BINNED_DATA_START"]["Fields"][1:])
            assert from_subject_samples == from_nmr_binned_data_samples

    if validate_factors:

        if "MS_METABOLITE_DATA" in mwtabfile:
            from_metabolite_data_factors = set(mwtabfile["MS_METABOLITE_DATA"]["MS_METABOLITE_DATA_START"]["Factors"])
            assert from_subject_factors == from_metabolite_data_factors


def validate_section(section, schema):
    """Validate section of ``mwTab`` formatted file.

    :param section: Section of :class:`~mwtab.mwtab.MWTabFile`.
    :param schema: Schema definition.
    :return: Validated section.
    :rtype: :py:class:`collections.OrderedDict`
    """
    return schema.validate(section)


def validate_file(mwtabfile, section_schema_mapping=section_schema_mapping, validate_samples=True, validate_factors=True):
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
    validated_mwtabfile = OrderedDict()

    try:
        _validate_samples_factors(mwtabfile, validate_samples, validate_factors)
    except Exception:
        raise

    for section_key, section in mwtabfile.items():
        try:
            schema = section_schema_mapping[section_key]
            section = validate_section(section=section, schema=schema)
            validated_mwtabfile[section_key] = section
        except Exception:
            raise

    return validated_mwtabfile
