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
    mwtab download <input-value> [--to-path=<path>] [--context=<context>] [--input-item=<item>] [--output-item=<item>] [--output-format=<format>] [--validate] [--verbose]
    mwtab extract metadata <from-path> <output-path> <key> ... [--extraction-format=<format>] [--no-header]
    mwtab extract metabolites <from-path> <output-path> (<key> <value>) ... [--extraction-format=<format>] [--no-header]


Options:
    -h, --help                      Show this screen.
    --version                       Show version.
    --verbose                       Print what files are processing.
    --validate                      Validate the mwTab file.
    --from-format=<format>          Input file format, available formats: mwtab, json [default: mwtab].
    --to-format=<format>            Output file format, available formats: mwtab, json [default: json].
    --mw-rest=<url>                 URL to MW REST interface
                                    [default: https://www.metabolomicsworkbench.org/rest/study/analysis_id/{}/mwtab/txt].
    --context=<context>             Type of resource to access from MW REST interface, available contexts: study,
                                    compound, refmet, gene, protein, moverz, exactmass [default: study].
    --input-item=<item>
    --output-item=<item>
    --output-format=<format>        Format for item to be retrieved in, available formats: mwtab, json, etc.
    --extraction-format=<format>    File format for extracted data/metadata to be save in, available formats: csv, json
                                    [default: csv].
    --no-header                     Include header at teh top of csv formatted files.

    <output-path> can take a "-" which will use stdout.
"""

from . import fileio
from . import mwrest
from . import mwextract
from .converter import Converter
from .validator import validate_file
from .mwschema import section_schema_mapping

from os import getcwd
from os.path import join

import json


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

    elif cmdargs["extract"]:
        mwfile_generator = fileio.read_files(cmdargs["<from-path>"])
        if cmdargs["metabolites"]:
            metabolites_dict = mwextract.extract_metabolites(mwfile_generator, cmdargs)
            if cmdargs["<output-path>"] != "-":
                if cmdargs["--extraction-format"] == "csv":
                    mwextract.write_metabolites_csv(cmdargs["<output-path>"], metabolites_dict)
                else:
                    mwextract.write_json(cmdargs["<output-path>"], metabolites_dict)
            else:
                print(json.dumps(metabolites_dict, indent=4, cls=mwextract.SetEncoder))

        elif cmdargs["metadata"]:
            metadata = dict()
            for mwtabfile in mwfile_generator:
                extracted_values = mwextract.extract_metadata(mwtabfile, cmdargs)
                [metadata.setdefault(key, set()).update(val) for (key, val) in extracted_values.items()]
            if cmdargs["<output-path>"] != "-":
                if cmdargs["--extraction-format"] == "csv":
                    mwextract.write_metadata_csv(cmdargs["<output-path>"], metadata)
                else:
                    mwextract.write_json(cmdargs["<output-path>"], metadata)
            else:
                print(metadata)
