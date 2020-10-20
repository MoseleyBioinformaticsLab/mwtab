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
    mwtab download url <url> [--to-path=<path>] [--verbose]
    mwtab download study all [--to-path=<path>] [--input-item=<item>] [--output-format=<format>] [--mw-rest=<url>] [--validate] [--verbose]
    mwtab download study <input-value> [--to-path=<path>] [--input-item=<item>] [--output-item=<item>] [--output-format=<format>] [--mw-rest=<url>] [--validate] [--verbose]
    mwtab download (study | compound | refmet | gene | protein) <input-item> <input-value> <output-item> [--output-format=<format>] [--to-path=<path>] [--mw-rest=<url>] [--verbose]
    mwtab download moverz <input-item> <m/z-value> <ion-type-value> <m/z-tolerance-value> [--to-path=<path>] [--mw-rest=<url>] [--verbose]
    mwtab download exactmass <LIPID-abbreviation> <ion-type-value> [--to-path=<path>] [--mw-rest=<url>] [--verbose]
    mwtab extract metadata <from-path> <to-path> <key> ... [--to-format=<format>] [--no-header]
    mwtab extract metabolites <from-path> <to-path> (<key> <value>) ... [--to-format=<format>] [--no-header]

Options:
    -h, --help                      Show this screen.
    --version                       Show version.
    --verbose                       Print what files are processing.
    --validate                      Validate the mwTab file.
    --from-format=<format>          Input file format, available formats: mwtab, json [default: mwtab].
    --to-format=<format>            Output file format [default: json].
                                    Available formats for convert:
                                        mwtab, json.
                                    Available formats for extract:
                                        json, csv.
    --mw-rest=<url>                 URL to MW REST interface
                                    [default: https://www.metabolomicsworkbench.org/rest/].
    --context=<context>             Type of resource to access from MW REST interface, available contexts: study,
                                    compound, refmet, gene, protein, moverz, exactmass [default: study].
    --input-item=<item>             Item to search Metabolomics Workbench with.
    --output-item=<item>            Item to be retrieved from Metabolomics Workbench.
    --output-format=<format>        Format for item to be retrieved in, available formats: mwtab, json [default: json]
    --no-header                     Include header at the top of csv formatted files.

    For extraction <to-path> can take a "-" which will use stdout.
"""

from . import fileio
from . import mwextract
from . import mwrest
from .converter import Converter
from .validator import validate_file
from .mwschema import section_schema_mapping

from os import getcwd, makedirs, path
from os.path import join
from urllib.parse import quote_plus

import json
import re


def check_filepath(filepath):
    """
    Method for validating that a given path directory exits. If not, the directory is created.

    :param str filepath: File path string.
    :return: None
    :rtype: :py:obj:`None`
    """
    if not path.exists(path.dirname(filepath)):
        dirname = path.dirname(filepath)
        if dirname:
            makedirs(dirname)


def download(context, cmdparams):
    """
    Method for creating Metabolomics Workbench REST URLs and requesting files based on given commandline arguments.
    Retrieved data is then saved out as specified.

    :param str context: String indicating the type of data to be accessed from the Metabolomics Workbench.
    :param dict cmdparams: Commandline arguments specifiying data to be accessed from Metabolomics Workbench.
    :return: None
    :rtype: :py:obj:`None`
    """
    mwresturl = mwrest.GenericMWURL({
        "context": context,
        "input_item": cmdparams["<input-item>"],
        "input_value": cmdparams["<input-value>"],
        "output_item": cmdparams["<output-item>"],
        "output_format": cmdparams["--output-format"],
    }).url
    mwrestfile = next(fileio.read_mwrest(mwresturl))

    # filepath validation
    extention = "txt"
    filepath = join(getcwd(), quote_plus(mwrestfile.source).replace(".", "_") + "." + (
        cmdparams["--output-format"] or extention))
    if cmdparams["--to-path"]:
        check_filepath(cmdparams["--to-path"])
        filepath = cmdparams["--to-path"]

    with open(filepath, "w", encoding="utf-8") as fh:
        mwrestfile.write(fh)


def cli(cmdargs):

    fileio.VERBOSE = cmdargs["--verbose"]
    mwrest.VERBOSE = cmdargs["--verbose"]
    fileio.MWREST = cmdargs["--mw-rest"]

    # mwtab convert ...
    if cmdargs["convert"]:
        converter = Converter(from_path=cmdargs["<from-path>"],
                              to_path=cmdargs["<to-path>"],
                              from_format=cmdargs["--from-format"],
                              to_format=cmdargs["--to-format"],
                              validate=cmdargs["--validate"])
        converter.convert()

    # mwtab validate ...
    elif cmdargs["validate"]:
        for mwfile in fileio.read_files(cmdargs["<from-path>"], validate=cmdargs["--validate"]):
            validate_file(mwtabfile=mwfile,
                          section_schema_mapping=section_schema_mapping,
                          validate_samples=True,
                          validate_factors=True,
                          verbose=True)

    # mwtab download ...
    elif cmdargs["download"]:

        # mwtab download url ...
        if cmdargs["<url>"]:
            mwrestfile = next(fileio.read_mwrest(cmdargs["<url>"]))
            with open(
                    cmdargs["--to-path"] or join(getcwd(), quote_plus(cmdargs["<url>"]).replace(".", "_") + "." +
                                                           cmdargs["--output-format"]),
                    "w"
            ) as fh:
                mwrestfile.write(fh)

        # mwtab download study ...
        elif cmdargs["study"]:

            # mwtab download study all ...
            if cmdargs["all"]:
                # mwtab download study all --input-item=analysis_id
                if not cmdargs["--input-item"] or cmdargs["--input-item"] == "analysis_id":
                    mwtabfiles = fileio.read_files(
                        *mwrest.generate_mwtab_urls(mwrest.analysis_ids(), cmdargs["--output-format"]),
                        valdate=cmdargs.get("--validate")
                    )
                    for mwtabfile in mwtabfiles:
                        print(quote_plus(mwtabfile.source).replace(".", "_"))
                        with open(cmdargs["--to-path"] or join(getcwd(), quote_plus(mwtabfile.source).replace(".", "_")) + "." + cmdargs["--output-format"], "w") as fh:
                            mwtabfile.write(fh, cmdargs["--output-format"])
                # mwtab download study all --input-item=study_id
                elif cmdargs["--input-item"] == "study_id":
                    mwtabfiles = fileio.read_files(
                        *mwrest.generate_mwtab_urls(mwrest.study_ids(), cmdargs["--output-format"]),
                        valdate=cmdargs.get("--validate")
                    )
                    for mwtabfile in mwtabfiles:
                        with open(cmdargs["--to-path"] or join(getcwd(), quote_plus(mwtabfile.source).replace(".", "_")) + "." + cmdargs["--output-format"], "w") as fh:
                            mwtabfile.write(fh, cmdargs["--output-format"])

            # mwtab download study <input_value> ...
            elif cmdargs["<input-value>"] and not cmdargs["<input-item>"]:
                input_item = cmdargs.get("--input-item")
                input_value = cmdargs["<input-value>"]
                if not input_item:
                    if input_value.isdigit():
                        input_value = "AN{}".format(input_value.zfill(6))
                        input_item = "analysis_id"
                    elif re.match(r'(AN[0-9]{6}$)', input_value):
                        input_item = "analysis_id"
                    elif re.match(r'(ST[0-9]{6}$)', input_value):
                        input_item = "study_id"
                mwresturl = mwrest.GenericMWURL({
                    "context": "study",
                    "input_item": input_item,
                    "input_value": input_value,
                    "output_item": cmdargs.get("--output-item") or "mwtab",
                    "output_format": cmdargs["--output-format"],
                }, cmdargs["--mw-rest"]).url
                mwrestfile = next(fileio.read_mwrest(mwresturl))
                with open(cmdargs["--to-path"] or join(getcwd(),
                                                       quote_plus(mwrestfile.source).replace(".", "_") + "." + cmdargs[
                                                           "--output-format"]),
                          "w") as fh:
                    mwrestfile.write(fh)

            # mwtab download (study | ...) <input_item> ...
            elif cmdargs["<input-item>"]:
                download("study", cmdargs)

        # mwtab download (... compound | refmet | gene | protein) ...
        elif cmdargs["compound"]:
            download("compound", cmdargs)
        elif cmdargs["refmet"]:
            download("refmet", cmdargs)
        elif cmdargs["gene"]:
            download("gene", cmdargs)
        elif cmdargs["protein"]:
            download("protein", cmdargs)

        # mwtab download moverz <input-value> <m/z-value> <ion-type-value> <m/z-tolerance-value> [--verbose]
        elif cmdargs["moverz"]:
            mwresturl = mwrest.GenericMWURL({
                "context": "moverz",
                "input_item": cmdargs["<input-item>"],
                "m/z_value": cmdargs["<m/z-value>"],
                "ion_type_value": cmdargs["<ion-type-value>"],
                "m/z_tolerance_value": cmdargs["<m/z-tolerance-value>"],
            }).url
            mwrestfile = next(fileio.read_mwrest(mwresturl))
            with open(cmdargs["--to-path"] or join(getcwd(), quote_plus(mwrestfile.source).replace(".", "_") + ".txt"),
                      "w") as fh:
                mwrestfile.write(fh)

        # mwtab download exactmass <LIPID-abbreviation> <ion-type-value> [--verbose]
        elif cmdargs["exactmass"]:
            mwresturl = mwrest.GenericMWURL({
                "context": "exactmass",
                "LIPID_abbreviation": cmdargs["<LIPID-abbreviation>"],
                "ion_type_value": cmdargs["<ion-type-value>"],
            }).url
            mwrestfile = next(fileio.read_mwrest(mwresturl))
            with open(cmdargs["--to-path"] or join(getcwd(), quote_plus(mwrestfile.source).replace(".", "_") + ".txt"),
                      "w") as fh:
                mwrestfile.write(fh)

    # mwtab extract ...
    elif cmdargs["extract"]:
        mwfile_generator = fileio.read_files(cmdargs["<from-path>"])
        if cmdargs["metabolites"]:
            metabolites_dict = mwextract.extract_metabolites(
                mwfile_generator,
                mwextract.generate_matchers([(
                    cmdargs["<key>"][i],
                    cmdargs["<value>"][i])
                    for i in range(len(cmdargs["<key>"]))])
            )

            if cmdargs["<to-path>"] != "-":
                if cmdargs["--to-format"] == "csv":
                    mwextract.write_metabolites_csv(cmdargs["<to-path>"], metabolites_dict, cmdargs["--no-header"])
                else:
                    mwextract.write_json(cmdargs["<to-path>"], metabolites_dict)
            else:
                print(json.dumps(metabolites_dict, indent=4, cls=mwextract.SetEncoder))

        elif cmdargs["metadata"]:
            metadata = dict()
            for mwtabfile in mwfile_generator:
                extracted_values = mwextract.extract_metadata(mwtabfile, cmdargs["<key>"])
                [metadata.setdefault(key, set()).update(val) for (key, val) in extracted_values.items()]
            if cmdargs["<to-path>"] != "-":
                if cmdargs["--to-format"] == "csv":
                    mwextract.write_metadata_csv(cmdargs["<to-path>"], metadata, cmdargs["--no-header"])
                else:
                    mwextract.write_json(cmdargs["<to-path>"], metadata)
            else:
                print(metadata)
