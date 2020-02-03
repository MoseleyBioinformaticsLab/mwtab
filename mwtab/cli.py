#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
The mwtab command-line interface
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Usage:
    mwtab -h | --help
    mwtab --version
    mwtab convert (<from-path> <to-path>) [--from-format=<format>] [--to-format=<format>] [--validate] [--mw-rest=<url>] [--verbose]
    mwtab validate <from-path> [--mw-rest=<url>] [--verbose]
    mwtab download <input-value> [--to-path=<path>] [--context=<context>] [--input-item=<item>] [--output-item=item(s)] [--output-format=<format>] [--validate] [--verbose]
    mwtab search-metabolites <from-path> <subject-type> <subject-species>

Options:
    -h, --help                      Show this screen.
    --version                       Show version.
    --verbose                       Print what files are processing.
    --validate                      Validate the mwTab file.
    --from-format=<format>          Input file format, available formats: mwtab, json [default: mwtab].
    --to-format=<format>            Output file format, available formats: mwtab, json [default: json].
    --mw-rest=<url>                 URL to MW REST interface [default:
                                    https://www.metabolomicsworkbench.org/rest/study/analysis_id/{}/mwtab/txt].
    --context=<context>             Type of resource to access from MW REST interface, available contexts: study,
                                    compound, refmet, gene, protein, moverz, exactmass [default: study].
    --input-item=<item>
    --output-item=<item(s)>
    --output-format=<format>        Format for item to be retrieved in, available formats: mwtab, json, etc.
"""

from . import fileio
from . import mwrest
from .converter import Converter
from .validator import validate_file
from .mwschema import section_schema_mapping

from os import getcwd
from os.path import join


def cli(cmdargs):

    fileio.VERBOSE = cmdargs["--verbose"]
    fileio.MWREST = cmdargs["--mw-rest"]

    if cmdargs["convert"]:
        converter = Converter(from_path=cmdargs["<from-path>"],
                              to_path=cmdargs["<to-path>"],
                              from_format=cmdargs["--from-format"],
                              to_format=cmdargs["--to-format"],
                              validate=cmdargs["--validate"])
        converter.convert()

    elif cmdargs["validate"]:
        for mwfile in fileio.read_files(cmdargs["<from-path>"], validate=cmdargs["--validate"]):
            validate_file(mwtabfile=mwfile,
                          section_schema_mapping=section_schema_mapping,
                          validate_samples=True,
                          validate_factors=True)

    elif cmdargs["download"]:
        if cmdargs["<input-value>"] == "all":
            an_ids = mwrest.analysis_ids()
            for mwfile in fileio.read_files(*an_ids):
                with open(join(cmdargs.get("--to-path") or getcwd(), mwfile.analysis_id+".txt"), "w") as outfile:
                    mwfile.write(outfile, "mwtab")
                    outfile.close()

        else:
            for mwfile in fileio.read_files(
                    mwrest.GenericMWURL(**{
                        "context": cmdargs.get("--context") or "study",
                        "input item": cmdargs.get("--input-item") or "analysis_id",
                        'input value': cmdargs["<input-value>"],
                        'output item': cmdargs.get("--output-item") or "mwtab",
                        'output format': cmdargs.get("--output-format") or "txt"
                    }).url):
                with open(join(cmdargs.get("--to-path") or getcwd(), mwfile.analysis_id+".txt"), "w") as outfile:
                    mwfile.write(outfile, "mwtab")
                    outfile.close()

    elif cmdargs["search-metabolites"]:
        metabolites = dict()
        analyses = set()

        for mwfile in fileio.read_files(cmdargs["<from-path>"]):
            if mwfile["SUBJECT"]["SUBJECT_TYPE"] == cmdargs["<subject-type>"] and \
                    mwfile["SUBJECT"]["SUBJECT_SPECIES"] == cmdargs["<subject-species>"]:
                if mwfile.get("METABOLITES"):
                    analyses.add(mwfile.analysis_id)
                    for metabolite in mwfile["METABOLITES"]["METABOLITES_START"]["DATA"]:
                        if metabolite["metabolite_name"] in metabolites.keys():
                            metabolites[metabolite["metabolite_name"]].append(mwfile.analysis_id)
                        else:
                            metabolites[metabolite["metabolite_name"]] = [mwfile.analysis_id]

        print("{} matched analyses:\n\t{}".format(len(analyses), analyses))

        metabolites = sorted([item for item in metabolites.items()], key=lambda x: len(x[1]), reverse=True)
        for metabolite in metabolites:
            print(metabolite[0])
            print("\t", metabolite[1])
