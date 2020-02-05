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


def extract_metabolites(mwtabfile, **kwargs):
    # metabolites = dict()
    # analyses = set()
    #
    # for mwfile in fileio.read_files(cmdargs["<from-path>"]):
    #     if mwfile["SUBJECT"]["SUBJECT_TYPE"] == cmdargs["<subject-type>"] and \
    #             mwfile["SUBJECT"]["SUBJECT_SPECIES"] == cmdargs["<subject-species>"]:
    #         if mwfile.get("METABOLITES"):
    #             analyses.add(mwfile.analysis_id)
    #             for metabolite in mwfile["METABOLITES"]["METABOLITES_START"]["DATA"]:
    #                 if metabolite["metabolite_name"] in metabolites.keys():
    #                     metabolites[metabolite["metabolite_name"]].append(mwfile.analysis_id)
    #                 else:
    #                     metabolites[metabolite["metabolite_name"]] = [mwfile.analysis_id]
    #
    # print("{} matched analyses:\n\t{}".format(len(analyses), analyses))
    #
    # metabolites = sorted([item for item in metabolites.items()], key=lambda x: len(x[1]), reverse=True)
    # for metabolite in metabolites:
    #     print(metabolite[0])
    #     print("\t", metabolite[1])
    pass


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


def write_json(output_path, extracted_dict):
    """
        For metabolites:
        metabolite_name, num_studies, num_analyses, num_samples,

        For metadata:
        key, value_1, ...
        :param str output_path:
        :param dict extracted_dict:
        :return:
        """
    with open(output_path+".json", "w") as outfile:
        for key in extracted_dict:
            extracted_dict[key] = list(extracted_dict[key])
        json.dump(extracted_dict, outfile, sort_keys=True, indent=4)
