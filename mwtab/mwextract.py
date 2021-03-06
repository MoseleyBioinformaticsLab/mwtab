#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
mwtab.mwextract
~~~~~~~~~~~

This module provides a number of functions and classes for extracting and saving data and metadata
stored in ``mwTab`` formatted files in the form of :class:`~mwtab.mwtab.MWTabFile`.
"""
import csv
import json
import os
import re


class ItemMatcher(object):
    """ItemMatcher class that can be called to match items from ``mwTab`` formatted files in the form of
    :class:`~mwtab.mwtab.MWTabFile`.
    """

    section_conversion = {
        "PR": "PROJECT",
        "ST": "STUDY",
        "SU": "SUBJECT",
        "CO": "COLLECTION",
        "TR": "TREATMENT",
        "SP": "SAMPLEPREP",
        "CH": "CHROMATOGRAPHY",
        "AN": "ANALYSIS",
        "MS": "MS",
        "NM": "NMR",
    }

    def __init__(self, full_key, value_comparison):
        """ItemMatcher initializer.

        :param str full_key: Key to match in :class:`~mwtab.mwtab.MWTabFile`.
        :param value_comparison: Value to match in :class:`~mwtab.mwtab.MWTabFile`.
        :type value_comparison: :class:`re.Pattern` or :py:class:`str`
        """
        self.full_key = full_key
        self.section, self.key = self.full_key.split(":")
        self.section = ItemMatcher.section_conversion[self.section]
        self.value_comparison = value_comparison

    def __call__(self, mwtabfile):
        """Match key value pair in :class:`~mwtab.mwtab.MWTabFile`.

        :param mwtabfile: Instance of :class:`~mwtab.mwtab.MWTabFile`.
        :type mwtabfile: :class:`~mwtab.mwtab.MWTabFile`
        :return: True if key and value are present, False otherwise.
        :rtype: :py:obj:`True` or :py:obj:`False`
        """
        return mwtabfile[self.section][self.key] == self.value_comparison


class ReGeXMatcher(ItemMatcher):
    """ReGeXMatcher class that can be called to match items from ``mwTab`` formatted files in the form of
    :class:`~mwtab.mwtab.MWTabFile` using regular expressions.
    """

    def __init__(self, full_key, value_comparison):
        """ItemMatcher initializer.

        :param str full_key: Key to match in :class:`~mwtab.mwtab.MWTabFile`.
        :param value_comparison: Value, in the form of a regular expression, to match in
        :class:`~mwtab.mwtab.MWTabFile`.
        :type value_comparison: :class:`re.Pattern`
        """
        super(ReGeXMatcher, self).__init__(full_key, value_comparison)

    def __call__(self, mwtabfile):
        """Match key value pair in :class:`~mwtab.mwtab.MWTabFile`.

        :param mwtabfile: Instance of :class:`~mwtab.mwtab.MWTabFile`.
        :type mwtabfile: :class:`~mwtab.mwtab.MWTabFile`
        :return: True if key and value are present, False otherwise.
        :rtype: :py:obj:`True` or :py:obj:`False`
        """
        return re.search(self.value_comparison, mwtabfile[self.section][self.key])


def generate_matchers(items):
    """Construct a generator that yields Matchers :class:`~mwtab.mwtab.ItemMatcher` or
    :class:`~mwtab.mwtab.ReGeXMatcher`.

    :param iterable items: Iterable object containing key value pairs to match.
    :return: Yields a Matcher object for each given item.
    :rtype: :class:`~mwtab.mwtab.ItemMatcher` or :class:`~mwtab.mwtab.ReGeXMatcher`
    """
    for item in items:
        if type(item[1]) == re.Pattern:
            yield ReGeXMatcher(item[0], item[1])

        else:
            yield ItemMatcher(item[0], item[1])


def extract_metabolites(sources, matchers):
    """Extract metabolite data from ``mwTab`` formatted files in the form of :class:`~mwtab.mwtab.MWTabFile`.

    :param generator sources: Generator of mwtab file objects (:class:`~mwtab.mwtab.MWTabFile`).
    :param generator matchers: Generator of matcher objects (:class:`~mwtab.mwextract.ItemMatcher` or
    :class:`~mwtab.mwextract.ReGeXMatcher`).
    :return: Extracted metabolites dictionary.
    :rtype: :py:class:`dict`
    """
    metabolites = dict()
    for mwtabfile in sources:
        if all(matcher(mwtabfile) for matcher in matchers):
            data_section_key = list(set(mwtabfile.keys()) & {"MS_METABOLITE_DATA", "NMR_METABOLITE_DATA", "NMR_BINNED_DATA"})[0]
            for data_list in mwtabfile[data_section_key]["Data"]:
                for test_key in (key for key in data_list.keys() if key != "Metabolite"):
                    try:
                        if float(data_list[test_key]) > 0:
                            metabolites.setdefault(data_list["Metabolite"], dict())\
                                .setdefault(mwtabfile.study_id, dict())\
                                .setdefault(mwtabfile.analysis_id, set())\
                                .add(test_key)
                    except Exception as e:
                        pass
    return metabolites


def extract_metadata(mwtabfile, keys):
    """Extract metadata data from ``mwTab`` formatted files in the form of :class:`~mwtab.mwtab.MWTabFile`.

    :param mwtabfile: mwTab file object for metadata to be extracted from.
    :type mwtabfile: :class:`~mwtab.mwtab.MWTabFile`
    :param list keys: List of metadata field keys for metadata values to be extracted.
    :return: Extracted metadata dictionary.
    :rtype: :py:class:`dict`
    """
    extracted_values = {}
    for section in mwtabfile:
        for metadata in mwtabfile[section]:
            for key in keys:
                if metadata == key:  # TODO: Allow for partial match, ReGeX, etc.
                    extracted_values.setdefault(key, set()).add(mwtabfile[section][metadata])

    return extracted_values


def write_metadata_csv(to_path, extracted_values, no_header=False):
    """Write extracted metadata :py:class:`dict` into csv file.

    Example:
    "metadata","value1","value2"
    "SUBJECT_TYPE","Human","Plant"

    :param str to_path: Path to output file.
    :param dict extracted_values: Metadata dictionary to be saved.
    :param bool no_header: If true header is not included, otherwise header is included.
    :return: None
    :rtype: :py:obj:`None`
    """
    if not os.path.exists(os.path.dirname(os.path.splitext(to_path)[0])):
        dirname = os.path.dirname(to_path)
        if dirname:
            os.makedirs(dirname)

    if not os.path.splitext(to_path)[1]:
        to_path += ".csv"

    with open(to_path, "w", newline="") as outfile:
        wr = csv.writer(outfile, quoting=csv.QUOTE_ALL)
        if not no_header:
            max_value_num = max([len(extracted_values[key]) for key in extracted_values.keys()])
            line_list = ["metadata"]
            line_list.extend(["value{}".format(num) for num in range(max_value_num)])
            wr.writerow(line_list)
        for key in extracted_values:
            line_list = [key]
            line_list.extend([val for val in sorted(extracted_values[key])])
            wr.writerow(line_list)


def write_metabolites_csv(to_path, extracted_values, no_header=False):
    """Write extracted metabolites data :py:class:`dict` into csv file.

    Example:
    "metabolite_name","num-studies","num_analyses","num_samples"
    "1,2,4-benzenetriol","1","1","24"
    "1-monostearin","1","1","24"
    ...

    :param str to_path: Path to output file.
    :param dict extracted_values: Metabolites data dictionary to be saved.
    :param bool no_header: If true header is not included, otherwise header is included.
    :return: None
    :rtype: :py:obj:`None`
    """
    csv_list = []
    for metabolite_key in extracted_values.keys():
        num_analyses = 0
        num_samples = 0
        for study_key in extracted_values[metabolite_key]:
            num_analyses += len(extracted_values[metabolite_key][study_key])
            for analysis_key in extracted_values[metabolite_key][study_key]:
                num_samples += len(extracted_values[metabolite_key][study_key][analysis_key])

        csv_list.append([
            metabolite_key,
            len(extracted_values[metabolite_key]),
            num_analyses,
            num_samples
        ])

    if not os.path.exists(os.path.dirname(os.path.splitext(to_path)[0])):
        dirname = os.path.dirname(to_path)
        if dirname:
            os.makedirs(dirname)

    if not os.path.splitext(to_path)[1]:
        to_path += ".csv"

    with open(to_path, "w", newline="") as outfile:
        wr = csv.writer(outfile, quoting=csv.QUOTE_ALL)
        if not no_header:
            wr.writerow(["metabolite_name", "num-studies", "num_analyses", "num_samples"])
        for line_list in csv_list:
            wr.writerow(line_list)


class SetEncoder(json.JSONEncoder):
    """SetEncoder class for encoding Python sets :py:class:`set` into json serializable objects :py:class:`list`.
    """

    def default(self, obj):
        """Method for encoding Python objects. If object passed is a set, converts the set to JSON serializable lists
        or calls base implementation.

        :param object obj: Python object to be json encoded.
        :return: JSON serializable object.
        :rtype: :py:class:`dict`, :py:class:`list`,
                :py:class:`tuple`, :py:class:`str`,
                :py:class:`int`, :py:class:`float`,
                :py:obj:`bool`, or :py:obj:`None`
        """
        if isinstance(obj, set):
            return list(obj)
        return json.JSONEncoder.default(self, obj)


def write_json(to_path, extracted_dict):
    """Write extracted data or metadata :py:class:`dict` into json file.

    Metabolites example:
    {
        "1,2,4-benzenetriol": {
            "ST000001": {
                "AN000001": [
                    "LabF_115816",
                    ...
                ]
            }
        }
    }

    Metadata example:
    {
        "SUBJECT_TYPE": [
            "Plant",
            "Human"
        ]
    }

    :param str to_path: Path to output file.
    :param dict extracted_dict: Metabolites data or metadata dictionary to be saved.
    :return: None
    :rtype: :py:obj:`None`
    """
    if not os.path.exists(os.path.dirname(os.path.splitext(to_path)[0])):
        dirname = os.path.dirname(to_path)
        if dirname:
            os.makedirs(dirname)

    if not os.path.splitext(to_path)[1]:
        to_path += ".json"

    with open(to_path, "w") as outfile:
        json.dump(extracted_dict, outfile, sort_keys=True, indent=4, cls=SetEncoder)
