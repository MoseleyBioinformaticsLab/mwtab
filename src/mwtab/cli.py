#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Note that this docstring has Sphinx directives in it. Those are removed in __main__ before passing to docopt.
"""
The mwtab Command Line Interface
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. code-block:: text

    Usage:
        mwtab -h | --help
        mwtab --version
        mwtab convert (<from-path> <to-path>) [--from-format=<format>] [--to-format=<format>] [--mw-rest=<url>] [--force] [--verbose]
        mwtab validate <from-path> [--to-path=<path>] [--mw-rest=<url>] [--force] [--silent]
        mwtab download url <url> [--to-path=<path>] [--verbose]
        mwtab download study all [--to-path=<path>] [--input-item=<item>] [--output-format=<format>] [--mw-rest=<url>] [--verbose]
        mwtab download study <input-value> [--to-path=<path>] [--input-item=<item>] [--output-item=<item>] [--output-format=<format>] [--mw-rest=<url>] [--verbose]
        mwtab download (study | compound | refmet | gene | protein) <input-item> <input-value> <output-item> [--output-format=<format>] [--to-path=<path>] [--mw-rest=<url>] [--verbose]
        mwtab download moverz <input-item> <m/z-value> <ion-type-value> <m/z-tolerance-value> [--to-path=<path>] [--mw-rest=<url>] [--verbose]
        mwtab download exactmass <LIPID-abbreviation> <ion-type-value> [--to-path=<path>] [--mw-rest=<url>] [--verbose]
        mwtab extract metadata <from-path> <to-path> <key> ... [--to-format=<format>] [--no-header] [--force]
        mwtab extract metabolites <from-path> <to-path> (<key> <value>) ... [--to-format=<format>] [--no-header] [--force]
    
    Options:
        -h, --help                           Show this screen.
        --version                            Show version.
        --verbose                            Print what files are processing.
        --silent                             Silence all standard output.
        --from-format=<format>               Input file format, available formats: mwtab, json [default: mwtab].
        --to-format=<format>                 Output file format [default: json].
                                             Available formats for convert:
                                                 mwtab, json.
                                             Available formats for extract:
                                                 json, csv.
        --mw-rest=<url>                      URL to MW REST interface
                                                [default: https://www.metabolomicsworkbench.org/rest/].
        --to-path=<path>                     Directory to save outputs into. Defaults to the current working directory.
                                             For the validate command, if the given path ends in '.json', then 
                                             all JSON file outputs will be condensed into that 1 file. Also for 
                                             the validate command no output files are saved unless this option is given.
        --prefix=<prefix>                    Prefix to add at the beginning of the output file name. Defaults to no prefix.
        --suffix=<suffix>                    Suffix to add at the end of the output file name. Defaults to no suffix.
        --context=<context>                  Type of resource to access from MW REST interface, available contexts: study,
                                             compound, refmet, gene, protein, moverz, exactmass [default: study].
        --input-item=<item>                  Item to search Metabolomics Workbench with.
        --output-item=<item>                 Item to be retrieved from Metabolomics Workbench.
        --output-format=<format>             Format for item to be retrieved in, available formats: mwtab, json.
        --no-header                          Include header at the top of csv formatted files.
        --force                              Ignore non-dictionary values in METABOLITES_DATA, METABOLITES, and EXTENDED tables for JSON files.
    
        For extraction <to-path> can take a "-" which will use stdout.
        All <from-path>'s can be single files, directories, or URLs.
    
    Documentation webpage: https://moseleybioinformaticslab.github.io/mwtab/
    GitHub webpage: https://github.com/MoseleyBioinformaticsLab/mwtab
"""

from os import getcwd
from os.path import join, isfile
from urllib.parse import quote_plus
import traceback
import json
import re
import sys
import time
import datetime
import pathlib

from . import fileio, mwextract, mwrest
from .converter import Converter
from .validator import validate_file
from .mwschema import ms_required_schema, nmr_required_schema
from .mwtab import MWTabFile


OUTPUT_FORMATS = {
    "txt": "txt",
    "mwtab": "txt",
    "json": "json",
    None: None
}
VERBOSE = False
# Note that 'url' is not a context for the Metabolomics Workbench REST API, but the code works better to have it in this list.
CONTEXTS = ['study', 'compound', 'refmet', 'gene', 'protein', 'moverz', 'exactmass', 'url']



def download(rest_params: dict, mwrest_base_url: str = mwrest.BASE_URL, full_url: str|None = None) -> mwrest.MWRESTFile:
    """Create Metabolomics Workbench REST URL and request file.
    
    Args:
        rets_params: A dictionary with values corresponding to the keywords in the Metabolomics Workbench REST specification.
                     For instance <context> would correspond to 'context' and <input item> would correspond to 'input_item'.
        mwrest_base_url: String for the base URL to use for accessing the Metabolomics Workbench REST interface.
        full_url: String representing a fully constructed URL to a Metabolomics Workbench REST endpoint. 
                  If given, all other parameters are ignored and this URL is used to download.
    
    Returns:
        The downloaded file as an MWRESTFile object if no exceptions occured, otherwise None with a raised exception.
    """
    if full_url:
        mwrestfile = next(fileio.read_mwrest(full_url))
    else:
        # create and validate a callable URL to pull data from Metabolomics Workbench's REST API
        fileurl = mwrest.GenericMWURL(rest_params, mwrest_base_url).url
        mwrestfile = next(fileio.read_mwrest(fileurl))
    
    return mwrestfile


def save_mwrest_file(mwrestfile: mwrest.MWRESTFile, to_path: str|None = None, output_format: str = 'txt') -> bool:
    """Save the given MWRESTFile object to the given path if it has text.
    
    Args:
        mwrestfile: The MWRESTFile object to save.
        to_path: The path to save the file to. Defaults to the current working directory.
        output_format: The format to save the file to. Should be 'txt', 'json', or 'mwtab'.
    
    Returns:
        True if the mwrestfile had text and was therefore saved, False otherwise.
    """
    if mwrestfile.text:  # if the text file isn't blank
        filename = quote_plus(mwrestfile.source).replace(".", "_")
        extension = OUTPUT_FORMATS[output_format]
        if to_path:
            fileio._create_save_path(to_path)
            if pathlib.Path(to_path).suffix:
                path_to_save = to_path
            else:
                path_to_save = join(to_path, filename + "." + extension)
        else:
            path_to_save = join(getcwd(), filename + "." + extension)
        
        with open(path_to_save, "w", encoding="utf-8") as fh:
            mwrestfile.write(fh)
        return True
    return False


def download_and_save_mwrest_file(rest_params: dict, to_path: str|None = None, 
                                  mwrest_base_url: str = mwrest.BASE_URL, full_url: str|None = None) -> None:
    """DRY function to combine downloading, saving, and error printing."""
    mwrestfile = download(rest_params, mwrest_base_url, full_url)
    extension = output_format if (output_format := rest_params.get('output_format')) else 'txt'
    if not save_mwrest_file(mwrestfile, to_path, extension):
        value = full_url if full_url else rest_params['input_value']
        print(f'When trying to download a file for the value, "{value}", '
              'a blank file or an error was returned, so no file was created for it.')


def download_and_save_ID_list(rest_params: dict, id_list: list[tuple[str, str]], verbose: bool,
                              to_path: str|None = None, 
                              mwrest_base_url: str = mwrest.BASE_URL, full_url: str|None = None) -> None:
    """Download and save a list of study and/or analysis IDs.
    
    Args:
        rest_params: A dictionary with values corresponding to the keywords in the Metabolomics Workbench REST specification.
                     For instance <context> would correspond to 'context' and <input item> would correspond to 'input_item'.
        id_list: A list of tuples that are pairs of IDs and their classification.
        verbose: If True, print more information about each ID being downloaded.
        to_path: The directory path to save the files to. Defaults to the current working directory.
        mwrest_base_url: String for the base URL to use for accessing the Metabolomics Workbench REST interface.
        full_url: String representing a fully constructed URL to a Metabolomics Workbench REST endpoint. 
                  If given, all other parameters are ignored and this URL is used to download.
    """
    for count, (input_id, input_item) in enumerate(id_list):
        rest_params['input_value'] = input_id
        rest_params['input_item'] = input_item
        if verbose:
            print("[{:4}/{:4}]".format(count+1, len(id_list)), input_id, datetime.datetime.now())
        try:
            download_and_save_mwrest_file(rest_params, to_path, mwrest_base_url, full_url)
        except Exception:
            print("Something went wrong and " + input_id + " could not be downloaded.")
            traceback.print_exc(file=sys.stdout)
            print()
        time.sleep(3)


def classify_input_value(input_value: str) -> tuple[str, str]:
    """Classify input_value as either 'analysis_id' or 'study_id'.
    
    If input_value is just a number, such as 000001, then 'AN' is added 
    to the front in the return value. The default classification is 'analysis_id'.
    
    Args:
        input_value: A string that should be a study ID or analysis ID.
    
    Returns:
        A tuple where the first element is the input_value, possibly modified, 
        and the second element is the classification, either 'study_id' or 'analysis_id'.
    """
    if input_value.isdigit():
        input_value = "AN{}".format(input_value.zfill(6))
        input_item = "analysis_id"
    elif re.match(r'(AN[0-9]{6}$)', input_value):
        input_item = "analysis_id"
    elif re.match(r'(ST[0-9]{6}$)', input_value):
        input_item = "study_id"
    else:
        input_item = "analysis_id"
    
    return input_value, input_item


def cli(cmdargs):
    """Implements the command line interface.

    param dict cmdargs: dictionary of command line arguments.
    """

    VERBOSE = cmdargs["--verbose"]
    force = cmdargs['--force']
    fileio.VERBOSE = cmdargs["--verbose"]
    silent = cmdargs['--silent']
    mwrest_base_url = cmdargs['--mw-rest']
    fileio.MWREST_URL = mwrest_base_url
    mwrest.BASE_URL = mwrest_base_url
    mwrest.VERBOSE = cmdargs["--verbose"]
    output_format = OUTPUT_FORMATS[cmdargs.get("--output-format")] if cmdargs.get("--output-format") else "txt"
    required_input_value = cmdargs.get('<input-value>')
    required_input_item = cmdargs.get('<input-item>')
    optional_input_item = cmdargs.get('--input-item')
    optional_to_path = cmdargs.get('--to-path')
    optional_output_item = cmdargs.get('--output-item')
    required_output_item = cmdargs.get('<output-item>')
    

    # mwtab convert ...
    if cmdargs["convert"]:
        converter = Converter(from_path=cmdargs["<from-path>"],
                              to_path=cmdargs["<to-path>"],
                              from_format=cmdargs["--from-format"],
                              to_format=cmdargs["--to-format"],
                              force=force)
        converter.convert()

    # mwtab validate ...
    elif cmdargs["validate"]:
        save_files = False
        consolidate_files = False
        consolidated_json = {}
        if optional_to_path:
            save_files = True
            fileio._create_save_path(optional_to_path)
            to_path = pathlib.Path(optional_to_path)
            if pathlib.Path(to_path).suffix == '.json':
                consolidate_files = True
        
        for i, (mwfile, e) in enumerate(fileio.read_with_class(cmdargs["<from-path>"], 
                                                               MWTabFile, 
                                                               {'duplicate_keys':True, 'force':force}, 
                                                               return_exceptions=True)):
            if e is not None:
                file_source = mwfile if isinstance(mwfile, str) else cmdargs["<from-path>"]
                print("Something went wrong when trying to read " + file_source)
                traceback.print_exception(e, file=sys.stdout)
                print()
                continue
            _, errors_list = validate_file(
                                            mwtabfile = mwfile,
                                            ms_schema = ms_required_schema, 
                                            nmr_schema = nmr_required_schema,
                                            verbose = not silent
                                            )
            if save_files:
                if consolidate_files:
                    consolidated_json[pathlib.Path(mwfile.source).stem] = errors_list
                else:
                    filename = pathlib.Path(mwfile.source).stem + '_validations.json'
                    with open(join(to_path, filename), "w", encoding="utf-8") as fh:
                        fh.write(json.dumps(errors_list, indent=2))
        
        if consolidate_files:
            with open(to_path, "w", encoding="utf-8") as fh:
                fh.write(json.dumps(consolidated_json, indent=2))
    
    # mwtab download ...
    elif cmdargs["download"]:
        context = [_context for _context in CONTEXTS if cmdargs.get(_context)][0]

        # mwtab download url ...
        if cmdargs["<url>"]:
            download_and_save_mwrest_file({}, optional_to_path, full_url = cmdargs['<url>'])

        # mwtab download study ...
        elif cmdargs["study"]:
            # mwtab download study all ...
            if cmdargs["all"]:
                # mwtab download study all ...
                # mwtab download study all --input-item=analysis_id ...
                # mwtab download study all --input-item=study_id ...
                if not optional_input_item or optional_input_item in ("analysis_id", "study_id"):
                    if not optional_input_item or optional_input_item == "analysis_id":
                        id_list = [(an_id, 'analysis_id') for an_id in mwrest.analysis_ids()]
                    elif optional_input_item == "study_id":
                        id_list = [(study_id, 'study_id') for study_id in mwrest.study_ids()]
                    rest_params = {'context': 'study',
                                   'output_item': 'mwtab',
                                   'output_format': output_format}
                    download_and_save_ID_list(rest_params, id_list, VERBOSE, optional_to_path, mwrest_base_url)

                else:
                    raise ValueError("Unknown \"--input-item\" {}".format(optional_input_item))

            # mwtab download study <input_value> ...
            elif required_input_value and not required_input_item:
                # Read a json file with a list and make a best attempt to parse and download studies and analyses from the list.
                if isfile(required_input_value):
                    with open(required_input_value, "r") as fh:
                        id_list = json.loads(fh.read())
                    
                    if len(id_list) > 1 and optional_to_path and pathlib.Path(optional_to_path).suffix:
                        print('Error: The given "--to-path" option is a path to a file, '
                              'but the given list of IDs to download is greater '
                              'than 1. Please specify a directory to save the files to.')
                        sys.exit(0)
                    
                    if optional_input_item:
                        if optional_input_item in ("analysis_id", "study_id"):
                            id_list = [(input_id, optional_input_item) for input_id in id_list]
                        else:
                            raise ValueError("Unknown \"--input-item\" {}".format(optional_input_item))
                    else:
                        id_list = [classify_input_value(_input_value) for _input_value in id_list]
                    
                    if VERBOSE:
                        print("Found {} Files to be Downloaded".format(len(id_list)))
                    
                    rest_params = {'context': 'study', 
                                    'output_item': optional_output_item if optional_output_item else 'mwtab',
                                    'output_format': output_format}
                    download_and_save_ID_list(rest_params, id_list, VERBOSE, optional_to_path, mwrest_base_url)

                # Assume input value is a single analysis or study id and use --input-item to decide which, default to analysis_id
                else:
                    input_item = cmdargs.get("--input-item")
                    input_value = required_input_value
                    if not input_item:
                        input_value, input_item = classify_input_value(input_value)
                    rest_params = {'context': 'study', 
                                   'input_value': input_value, 
                                   'input_item': input_item,
                                   'output_item': optional_output_item if optional_output_item else 'mwtab',
                                   'output_format': output_format}
                    download_and_save_mwrest_file(rest_params, optional_to_path, mwrest_base_url)

            # mwtab download (study | ...) <input_item> ...
            elif required_input_item:
                rest_params = {'context': 'study', 
                               'input_value': required_input_value, 
                               'input_item': required_input_item,
                               'output_item': required_output_item,
                               'output_format': output_format}
                download_and_save_mwrest_file(rest_params, optional_to_path, mwrest_base_url)
        
        # mwtab download (... | compound | refmet | gene | protein) ...
        elif context in ['compound', 'refmet', 'gene', 'protein']:
            rest_params = {'context': context, 
                           'input_value': required_input_value, 
                           'input_item': required_input_item,
                           'output_item': required_output_item,
                           'output_format': output_format}
            download_and_save_mwrest_file(rest_params, optional_to_path, mwrest_base_url)
                
        # mwtab download moverz <input-value> <m/z-value> <ion-type-value> <m/z-tolerance-value> [--verbose]
        elif cmdargs["moverz"]:
            rest_params = {'context': 'moverz', 
                           'input_item': required_input_item, 
                           "m/z_value": cmdargs["<m/z-value>"],
                           "ion_type_value": cmdargs["<ion-type-value>"],
                           "m/z_tolerance_value": cmdargs["<m/z-tolerance-value>"]}
            download_and_save_mwrest_file(rest_params, optional_to_path, mwrest_base_url)

        # mwtab download exactmass <LIPID-abbreviation> <ion-type-value> [--verbose]
        elif cmdargs["exactmass"]:
            rest_params = {'context': 'exactmass', 
                           "LIPID_abbreviation": cmdargs["<LIPID-abbreviation>"],
                           "ion_type_value": cmdargs["<ion-type-value>"],}
            download_and_save_mwrest_file(rest_params, optional_to_path, mwrest_base_url)

    # mwtab extract ...
    elif cmdargs["extract"]:
        mwfile_generator = fileio.read_with_class(cmdargs["<from-path>"], 
                                                  MWTabFile, 
                                                  {'duplicate_keys':True, "force": force}, 
                                                  return_exceptions=True)
        if cmdargs["metabolites"]:
            metabolites_dict = mwextract.extract_metabolites(
                mwfile_generator,
                mwextract.generate_matchers(
                    [(cmdargs["<key>"][i],
                      cmdargs["<value>"][i] if not cmdargs["<value>"][i][:2] == "r'" else re.compile(cmdargs["<value>"][i][2:-1]))
                     for i in range(len(cmdargs["<key>"]))]
                )
            )
            
            if metabolites_dict:
                if cmdargs["<to-path>"] != "-":
                    if cmdargs["--to-format"] == "csv":
                        mwextract.write_metabolites_csv(cmdargs["<to-path>"], metabolites_dict, cmdargs["--no-header"])
                    else:
                        mwextract.write_json(cmdargs["<to-path>"], metabolites_dict)
                else:
                    print(json.dumps(metabolites_dict, indent=4, cls=mwextract.SetEncoder))
            else:
                print("No metabolites extracted. No file was saved. " 
                      "This is likely due to key value pairs filtering all of the studies out. "
                      "Check your key value pairs and data.")

        elif cmdargs["metadata"]:
            metadata = dict()
            for mwtabfile, e in mwfile_generator:
                if e is not None:
                    file_source = mwtabfile if isinstance(mwtabfile, str) else cmdargs["<from-path>"]
                    print("Something went wrong when trying to read " + file_source)
                    traceback.print_exception(e, file=sys.stdout)
                    print()
                    continue
                extracted_values = mwextract.extract_metadata(mwtabfile, cmdargs["<key>"])
                [metadata.setdefault(key, set()).update(val) for (key, val) in extracted_values.items()]
            if metadata:
                if cmdargs["<to-path>"] != "-":
                    if cmdargs["--to-format"] == "csv":
                        mwextract.write_metadata_csv(cmdargs["<to-path>"], metadata, cmdargs["--no-header"])
                    else:
                        mwextract.write_json(cmdargs["<to-path>"], metadata)
                else:
                    print(metadata)
            else:
                print("No metadata extracted. No file was saved.")

