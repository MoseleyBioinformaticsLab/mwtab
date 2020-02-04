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
    mwtab extract metabolites <from-path> <output-filename> <key>=<value> ... [--extraction-format=<format>]
    mwtab extract metadata <from-path> <output-filename> <key> ... [--extraction-format=<format>]


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

    <output-filename> can take a "-" which will use stdout.
"""

from . import fileio
from . import mwrest
from . import mwextract
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

    elif cmdargs["extract"]:
        if cmdargs["metabolites"]:
            mwextract.extract_metabolites(cmdargs)
        elif cmdargs["metadata"]:
            mwextract.extract_metabolites(cmdargs)
