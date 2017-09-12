#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
The mwtab command-line interface
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Usage:
    mwtab -h | --help
    mwtab --version
    mwtab convert (<from-path> <to-path>) [--from-format=<format>] [--to-format=<format>] [--validate] [--verbose]
    mwtab validate <from-path> [--verbose]
    
Options:
    -h, --help                      Show this screen.
    --version                       Show version.
    --verbose                       Print what files are processing.
    --validate                      Validate the mwTab file.
    --from-format=<format>          Input file format, available formats: mwtab, json [default: mwtab].
    --to-format=<format>            Output file format, available formats: mwtab, json [default: json].
    --mw-rest=<url>                 URL to MW REST interface [default: "http://www.metabolomicsworkbench.org/rest/study/analysis_id/{}/mwtab/txt"].
"""

import mwtab
from mwtab.converter import Converter
from mwtab.validator import validate_file
from mwtab.mwschema import section_schema_mapping


def cli(cmdargs):

    mwtab.fileio.VERBOSE = cmdargs["--verbose"]
    mwtab.fileio.MWREST = cmdargs["--verbose"]

    if cmdargs["convert"]:
        converter = Converter(from_path=cmdargs["<from-path>"],
                              to_path=cmdargs["<to-path>"],
                              from_format=cmdargs["--from-format"],
                              to_format=cmdargs["--to-format"])
        converter.convert()

    elif cmdargs["validate"]:
        for mwfile in mwtab.read_files(cmdargs["<from-path>"]):
            validate_file(mwtabfile=mwfile,
                          section_schema_mapping=section_schema_mapping,
                          validate_samples=True,
                          validate_factors=True)
