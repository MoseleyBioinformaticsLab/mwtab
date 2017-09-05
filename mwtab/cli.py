#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
The mwtab command-line interface
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Usage:
    mwtab -h | --help
    mwtab --version
    mwtab convert (<from_path> <to_path>) [--from_format=<format>] [--to_format=<format>] [--validate] [--verbose]
    mwtab validate <from_path> [--verbose]
    
Options:
    -h, --help                      Show this screen.
    --version                       Show version.
    --verbose                       Print what files are processing.
    --validate                      Validate the mwTab file.
    --from_format=<format>          Input file format, available formats: mwtab, json [default: mwtab].
    --to_format=<format>            Output file format, available formats: mwtab, json [default: json].
"""

import mwtab
from mwtab.converter import Converter
from mwtab.validator import validate_file
from mwtab.mwschema import section_schema_mapping


def cli(cmdargs):

    if cmdargs["convert"]:
        converter = Converter(from_path=cmdargs["<from_path>"],
                              to_path=cmdargs["<to_path>"],
                              from_format=cmdargs["--from_format"],
                              to_format=cmdargs["--to_format"])
        converter.convert()

    elif cmdargs["validate"]:
        for mwfile in mwtab.read_files(cmdargs["<from_path>"]):
            validate_file(mwtabfile=mwfile,
                          section_schema_mapping=section_schema_mapping,
                          validate_samples=True,
                          validate_factors=True)
