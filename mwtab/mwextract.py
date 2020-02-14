#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
mwtab.mwextract
~~~~~~~~~~~

This module provides a number of functions for extracting metadata and data
stored in ``mwTab`` formatted files.
"""
import csv
import json
import re


class ItemMatcher(object):
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
        self.full_key = full_key
        self.section, self.key = self.full_key.split(":")
        self.section = ItemMatcher.section_conversion[self.section]
        self.value_comparison = value_comparison

    def __call__(self, mwtabfile):
        return mwtabfile[self.section][self.key] == self.value_comparison


class ReGeXMatcher(ItemMatcher):

    def __init__(self, full_key, value_comparison):
        super(ReGeXMatcher, self).__init__(full_key, value_comparison)

    def __call__(self, mwtabfile):
        return re.search(self.value_comparison, mwtabfile[self.section][self.key])


def extract_metabolites(mwfile_generator, kwargs):

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
    """

    :param mwtabfile:
    :param kwargs:
    :return:
    """
    extracted_values = {}
    for section in mwtabfile:
        for metadata in mwtabfile[section]:
            for key in kwargs["<key>"]:
                if metadata == key:  # TODO: Allow for partial match, ReGeX, etc.
                    extracted_values.setdefault(key, set()).add(mwtabfile[section][metadata])

    return extracted_values


def write_metadata_csv(output_path, extracted_values):
    """
    For metabolites:
    metabolite_name, num_studies, num_analyses, num_samples,

    For metadata:
    key, value_1, ...
    :param str output_path:
    :param dict extracted_values:
    :return:
    """
    with open(output_path+".csv", "w") as outfile:
        wr = csv.writer(outfile, quoting=csv.QUOTE_ALL)
        for key in extracted_values:
            line_list = [key]
            line_list.extend([val for val in sorted(extracted_values[key])])
            wr.writerow(line_list)


def write_metabolites_csv(output_path, extracted_values):
    """
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

    with open(output_path + ".csv", "w") as outfile:
        wr = csv.writer(outfile, quoting=csv.QUOTE_ALL)
        for line_list in csv_list:
            wr.writerow(line_list)


class SetEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, set):
            return list(obj)
        return json.JSONEncoder.default(self, obj)


def write_json(output_path, extracted_dict):
    """
    :param str output_path:
    :param dict extracted_dict:
    :return:
    """
    with open(output_path+".json", "w") as outfile:
        json.dump(extracted_dict, outfile, sort_keys=True, indent=4, cls=SetEncoder)
