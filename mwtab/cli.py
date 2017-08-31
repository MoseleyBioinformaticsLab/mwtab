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

def cli(cmdargs):

    if cmdargs["convert"]:
        pass