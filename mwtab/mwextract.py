#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
mwtab.mwextract
~~~~~~~~~~~

This module provides a number of functions for extracting metadata and data
stored in ``mwTab`` formatted files.
"""


def extract_metabolites(**kwargs):
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


def extract_metadata(kwargs):
    pass


def write_csv():
    pass


def write_json():
    pass
