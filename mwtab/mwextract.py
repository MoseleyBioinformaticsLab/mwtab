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

        :param str full_key:
        :param str value_comparison:
        """
        self.full_key = full_key
        self.section, self.key = self.full_key.split(":")
        self.section = ItemMatcher.section_conversion[self.section]
        self.value_comparison = value_comparison

    def __call__(self, mwtabfile):
        """Match key value pair in :class:`~mwtab.mwtab.MWTabFile`.

        :param mwtabfile:
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

        :param str full_key:
        :param str value_comparison:
        """
        super(ReGeXMatcher, self).__init__(full_key, value_comparison)

    def __call__(self, mwtabfile):
        """Match key value pair in :class:`~mwtab.mwtab.MWTabFile`.

        :param mwtabfile:
        :type mwtabfile: :class:`~mwtab.mwtab.MWTabFile`
        :return: True if key and value are present, False otherwise.
        :rtype: :py:obj:`True` or :py:obj:`False`
        """
        return re.search(self.value_comparison, mwtabfile[self.section][self.key])


def extract_metabolites(mwfile_generator, kwargs):
    """Extract metabolite data from ``mwTab`` formatted files in the form of :class:`~mwtab.mwtab.MWTabFile`.

    :param generator mwfile_generator:
    :param kwargs:
    :return: Extracted metabolites dictionary.
    :rtype: :py:class:`dict`
    """
    metabolites = dict()
    matchers = [ItemMatcher(kwargs["<key>"][i], kwargs["<value>"][i]) for i in range(len(kwargs["<key>"]))]
    for mwtabfile in mwfile_generator:
        if all(matcher(mwtabfile) for matcher in matchers):
            for metabolite in mwtabfile["METABOLITES"]["METABOLITES_START"]["DATA"]:
                for data_list in mwtabfile["MS_METABOLITE_DATA"]["MS_METABOLITE_DATA_START"]["DATA"]:
                    sample_keys = [k for k in data_list.keys() if k != "metabolite_name"]
                    for k in sample_keys:
                        if float(data_list[k]) > 0:
                            metabolites.setdefault(metabolite["metabolite_name"], dict())\
                                .setdefault(mwtabfile.study_id, dict())\
                                .setdefault(mwtabfile.analysis_id, set())\
                                .add(k)
    return metabolites


def extract_metadata(mwtabfile, kwargs):
    """Extract metadata data from ``mwTab`` formatted files in the form of :class:`~mwtab.mwtab.MWTabFile`.

    :param generator mwtabfile:
    :param kwargs:
    :return: Extracted metadata dictionary.
    :rtype: :py:class:`dict`
    """
    extracted_values = {}
    for section in mwtabfile:
        for metadata in mwtabfile[section]:
            for key in kwargs["<key>"]:
                if metadata == key:  # TODO: Allow for partial match, ReGeX, etc.
                    extracted_values.setdefault(key, set()).add(mwtabfile[section][metadata])

    return extracted_values


def write_metadata_csv(output_path, extracted_values, no_header=False):
    """Write extracted metadata :py:class:`dict` into csv file.

    Example:
    "metadata","value1","value2"
    "SUBJECT_TYPE","Human","Plant"

    :param str output_path:
    :param dict extracted_values:
    :param bool no_header:
    :return: None
    :rtype: :py:obj:`None`
    """
    with open(output_path+".csv", "w") as outfile:
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


def write_metabolites_csv(output_path, extracted_values, no_header=False):
    """Write extracted metabolites data :py:class:`dict` into csv file.


    :param str output_path:
    :param dict extracted_values:
    :param bool no_header:
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
    print(csv_list)

    with open(output_path + ".csv", "w") as outfile:
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


def write_json(output_path, extracted_dict):
    """Write extracted data or metadata :py:class:`dict` into json file.

    :param str output_path:
    :param dict extracted_dict:
    :return: None
    :rtype: :py:obj:`None`
    """
    with open(output_path+".json", "w") as outfile:
        json.dump(extracted_dict, outfile, sort_keys=True, indent=4, cls=SetEncoder)
