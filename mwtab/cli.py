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
    mwtab request <input-value> [--to-path=<path>] [--context=<context>] [--input-item=<item>] [--output-item=item(s)] [--output-format=<format>] [--validate] [--verbose]
    
Options:
    -h, --help                      Show this screen.
    --version                       Show version.
    --verbose                       Print what files are processing.
    --validate                      Validate the mwTab file.
    --from-format=<format>          Input file format, available formats: mwtab, json [default: mwtab].
    --to-format=<format>            Output file format, available formats: mwtab, json [default: json].
    --mw-rest=<url>                 URL to MW REST interface [default: https://www.metabolomicsworkbench.org/rest/study/analysis_id/{}/mwtab/txt].

    --to-path                       File path to directory for file to be saved in [default: cwd].
    --context=<context>             Type of resource to access from MW REST interface, available contexts: study, compound, refmet, gene, protein, moverz, exactmass [default: study].
    --input-item=<item>
    --output-item=<item(s)>
    --output-format=<format>        Format for item to be retrieved in, available formats: mwtab, json, etc. [default: mwtab].
"""

from . import fileio
from .converter import Converter
from .validator import validate_file
from .mwschema import section_schema_mapping


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

    elif cmdargs["request"]:

