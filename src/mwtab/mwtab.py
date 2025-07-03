#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
mwtab.mwtab
~~~~~~~~~~~

This module provides the :class:`~mwtab.mwtab.MWTabFile` class
that stores the data from a single ``mwTab`` formatted file in the
form of an :py:class:`~collections.OrderedDict`. Data can be accessed
directly from the :class:`~mwtab.mwtab.MWTabFile` instance using
bracket accessors.

The data is divided into a series of "sections" which each contain a
number of "key-value"-like pairs. Also, the file contains a specially
formatted ``SUBJECT_SAMPLE_FACTOR`` block and blocks of data between 
``*_START`` and ``*_END``.
"""

from __future__ import print_function, division, unicode_literals
from collections import OrderedDict
import io
import sys
import json
import re
import copy
from itertools import zip_longest

import json_duplicate_keys as jdks
import pandas

from .tokenizer import tokenizer
from .validator import validate_file
from .mwschema import section_schema_mapping
from .duplicates_dict import DuplicatesDict

DUPLICATE_KEY_REGEX = r'(.*)\{\{\{.*\}\}\}'

# The stuff before the MWTabFile class is all to do with being able to handle duplicate keys from a JSON file.
# Python's parser can't do it and you have to do some workarounds for it.

# From https://stackoverflow.com/questions/14902299/json-loads-allows-duplicate-keys-in-a-dictionary-overwriting-the-first-value
def _handle_duplicate_keys(ordered_pairs):
    """Use special type to store duplicate keys."""
    d = DuplicatesDict()
    for k, v in ordered_pairs:
        if k in d:
            d.set(k, v, ordered_dict = True)
        #     if isinstance(d[k], _duplicate_key_list):
        #         d[k].append(v)
        #     else:
        #         d[k] = _duplicate_key_list([d[k],v])
        # else:
        #    d[k] = v
    
    if len({k for k, v in ordered_pairs}) == len(ordered_pairs):
        return d.getObject()
    else:
        return d


SORT_KEYS = False
INDENT = 4
def _JSON_serializer_for_dupe_class(o):
    """ """
    return o.dumps(sort_keys=SORT_KEYS, indent=INDENT)


def _match_process(matchobj):
    temp_string = matchobj.group(2)
    temp_string = temp_string.replace('\\"', '"')
    temp_string = temp_string.replace('\\n',  '\n' + ' '*(3*INDENT))
    return matchobj.group(1) + ': {' + temp_string + '}'


def _parse_header_input(input_str):
    """Attempt to parse a header string into a dict.
    """
    match_re = re.match(r"#METABOLOMICS WORKBENCH( )?([^: ]+ )?([A-Z_]+:\w+ ?)*", input_str)
    if match_re:
        key_values = re.findall(r" (\w*:\w*)", input_str)
        temp_dict = {}
        for key_value in key_values:
            key, new_value = key_value.split(":")
            temp_dict[key] = new_value
        
        return temp_dict
    else:
        raise ValueError("header cannot be set because it is not of the form \"#METABOLOMICS WORKBENCH( )?([^: ]+ )?([A-Z_]+:\w+ ?)*\"")


# Descriptor to handle the convenience properties for MWTabFile.
# https://realpython.com/python-descriptors/
class MWTabProperty:
    def __set_name__(self, owner, name):
        self._name = name
    
    def __get__(self, obj, type=None):
        if self._name == "study_id" or self._name == "analysis_id":
            try:
                return obj["METABOLOMICS WORKBENCH"].get(self._name.upper())
            except Exception:
                return None
        
        if self._name == "header":
            try:
                header_str = "#METABOLOMICS WORKBENCH"
                pairs = " ".join([item[0] + ":" + item[1] for item in obj["METABOLOMICS WORKBENCH"].items() 
                                  if item[0] not in ["VERSION", "CREATED_ON"]])
                header_str += " " + pairs
                return header_str
            except Exception:
                return None
    
    def __set__(self, obj, value):
        if not isinstance(value, str):
            raise TypeError("The value for " + self._name + " must be a string.")
        
        if "METABOLOMICS WORKBENCH" in obj:
            if self._name == "study_id" or self._name == "analysis_id":
                if isinstance(obj["METABOLOMICS WORKBENCH"], (dict, OrderedDict)):
                    obj["METABOLOMICS WORKBENCH"][self._name.upper()] = value
                else:
                    raise TypeError("The \"METABOLOMICS WORKBENCH\" key is not a dictionary, so " + self._name + " cannot be set.")
            
            if self._name == "header":
                temp_dict = _parse_header_input(value)
                obj["METABOLOMICS WORKBENCH"].update(temp_dict)
        else:
            if self._name == "study_id" or self._name == "analysis_id":
                obj["METABOLOMICS WORKBENCH"] = {self._name.upper() : value}
            
            if self._name == "header":
                temp_dict = _parse_header_input(value)
                obj["METABOLOMICS WORKBENCH"]  = temp_dict
            
    def __delete__(self, obj):
        del obj.__dict__[self._name]


class MWTabFile(OrderedDict):
    """MWTabFile class that stores data from a single ``mwTab`` formatted file in
    the form of :py:class:`collections.OrderedDict`.
    """

    prefixes = {
        "METABOLOMICS WORKBENCH": "",
        "PROJECT": "PR:",
        "STUDY": "ST:",
        "SUBJECT": "SU:",
        "SUBJECT_SAMPLE_FACTORS": "",
        "COLLECTION": "CO:",
        "TREATMENT": "TR:",
        "SAMPLEPREP": "SP:",
        "CHROMATOGRAPHY": "CH:",
        "ANALYSIS": "AN:",
        "MS": "MS:",
        "NMR": "NM:",
        "NM": "NM:",
        "MS_METABOLITE_DATA": "",
        "NMR_METABOLITE_DATA": "",
        "NMR_BINNED_DATA": "",
        "METABOLITES": "",
    }
    
    result_file_keys = ["UNITS", "Has m/z", "Has RT", "RT units"]
    data_section_keys = {"MS_METABOLITE_DATA", "NMR_METABOLITE_DATA", "NMR_BINNED_DATA"}
    
    study_id = MWTabProperty()
    analysis_id = MWTabProperty()
    header = MWTabProperty()

    def __init__(self, source, compatability_mode=False, *args, **kwds):
        """File initializer.

        :param str source: Source a `MWTabFile` instance was created from.
        """
        super(MWTabFile, self).__init__(*args, **kwds)
        self.source = source
        self._factors = None
        self._samples = None
        self._metabolite_header = None
        self._extended_metabolite_header = None
        self._binned_header = None
        self._short_headers = set()
        self._duplicate_sub_sections = {}
        self.compatability_mode = compatability_mode
        if compatability_mode:
            self._default_dict_type = DuplicatesDict
        else:
            self._default_dict_type = OrderedDict
        # self._study_id = None
        # self._analysis_id = None
        # self._header = None
    
    @property
    def data_section_key(self):
        """Easily determine the data_section_key.
        
        The key will be one of "MS_METABOLITE_DATA", "NMR_METABOLITE_DATA", or "NMR_BINNED_DATA", 
        but will be None if none of those keys are found.
        """
        data_section_key = list(set(self.keys()) & self.data_section_keys)
        if data_section_key:
            return data_section_key[0]
        return None
    
    def set_table_as_pandas(self, df, table_name, clear_header=False):
        """Return the given table_name as a pandas.DataFrame.
        
        table_name must be one of "Metabolites", "Extended", or "Data".
        
        :param pandas.DataFrame df: pandas.DataFrame that will be used to update the object.
        :param str table_name: the name of the table to set from df.
        :param bool clear_header: if True, sets the appropriate header to None, otherwise sets it to the df columns.
        :return: None
        :rtype: :py:obj:`None`
        """
        data_section_key = self.data_section_key
        self[data_section_key][table_name] = [self._default_dict_type(data_dict) for data_dict in df.fillna('').astype(str).to_dict(orient='records')]
        if clear_header:
            if table_name == 'Metabolites':
                self._metabolite_header = None
            elif table_name == 'Extended':
                self._extended_metabolite_header = None
            elif table_name == 'Data':
                self._samples = None
                if "BINNED" in data_section_key:
                    self._binned_header = self._samples
        elif df.shape[1] > 0:
            if table_name == 'Metabolites':
                self._metabolite_header = [k for k in self[data_section_key][table_name][0].keys()][1:]
            elif table_name == 'Extended':
                self._extended_metabolite_header = [k for k in self[data_section_key][table_name][0].keys()][1:]
            elif table_name == 'Data':
                self._samples = [k for k in self[data_section_key][table_name][0].keys()][1:]
                if "BINNED" in data_section_key:
                    self._binned_header = self._samples
    
    def set_metabolites_from_pandas(self, df, clear_header=False):
        """Update MWTabFile based on provided pandas.DataFrame.
        
        Overwrite the current list of dicts in self[data_section_key]['Metabolites'] with the 
        values in df. Also overwrites self._metabolite_header with the columns in df, excluding the 
        first column. df is assumed to have the first column as the 'Metabolites' column.
        
        :param pandas.DataFrame df: pandas.DataFrame that will be used to update the object.
        :param bool clear_header: if True, sets _metabolite_header to None, otherwise sets it to the df columns.
        :return: None
        :rtype: :py:obj:`None`
        """
        self.set_table_as_pandas(df, 'Metabolites', clear_header)
        # data_section_key = self.data_section_key
        # self[data_section_key]['Metabolites'] = [self._default_dict_type(data_dict) for data_dict in df.fillna('').astype(str).to_dict(orient='records')]
        # if df.shape[1] > 0:
        #     self._metabolite_header = [k for k in self[data_section_key]['Metabolites'][0].keys()][1:]
        # elif clear_header:
        #     self._metabolite_header = None
    
    def set_extended_from_pandas(self, df, clear_header=False):
        """Update MWTabFile based on provided pandas.DataFrame.
        
        Overwrite the current list of dicts in self[data_section_key]['Extended'] with the 
        values in df. Also overwrites self._extended_metabolite_header with the columns in df, excluding the 
        first column. df is assumed to have the first column as the 'Metabolites' column.
        
        :param pandas.DataFrame df: pandas.DataFrame that will be used to update the object.
        :param bool clear_header: if True, sets _extended_metabolite_header to None, otherwise sets it to the df columns.
        :return: None
        :rtype: :py:obj:`None`
        """
        self.set_table_as_pandas(df, 'Extended', clear_header)
        # data_section_key = self.data_section_key
        # self[data_section_key]['Extended'] = [self._default_dict_type(data_dict) for data_dict in df.fillna('').astype(str).to_dict(orient='records')]
        # if df.shape[1] > 0:
        #     self._extended_metabolite_header = [k for k in self[data_section_key]['Extended'][0].keys()][1:]
        # elif clear_header:
        #     self._extended_metabolite_header = None
    
    def set_metabolites_data_from_pandas(self, df, clear_header=False):
        """Update MWTabFile based on provided pandas.DataFrame.
        
        Overwrite the current list of dicts in self[data_section_key]['Data'] with the 
        values in df. Also overwrites self._samples with the columns in df, excluding the 
        first column. df is assumed to have the first column as the 'Metabolites' column.
        
        :param pandas.DataFrame df: pandas.DataFrame that will be used to update the object.
        :param bool clear_header: if True, sets _samples to None, otherwise sets it to the df columns.
        :return: None
        :rtype: :py:obj:`None`
        """
        self.set_table_as_pandas(df, 'Data', clear_header)
        # data_section_key = self.data_section_key
        # self[data_section_key]['Data'] = [self._default_dict_type(data_dict) for data_dict in df.fillna('').astype(str).to_dict(orient='records')]
        # if df.shape[1] > 0:
        #     self._samples = [k for k in self[data_section_key]['Data'][0].keys()][1:]
        # elif clear_header:
        #     self._samples = None
        # if "BINNED" in data_section_key:
        #     self._binned_header = self._samples
    
    def get_table_as_pandas(self, table_name):
        """Return the given table_name as a pandas.DataFrame.
        
        table_name must be one of "Metabolites", "Extended", or "Data". Note that 
        if there are duplicate column names, they will have a string appended to the 
        end of the name like {{{_\d+_}}}.
        
        :param str table_name: the name of the table to return as a pandas.DataFrame.
        :return: The list of dicts for the given table_name as a pandas.DataFrame.
        :rtype: pandas.DataFrame
        """
        data_section_key = self.data_section_key
        if data_section_key:
            if self.compatability_mode:
                temp_list = [duplicates_dict._JSON_DUPLICATE_KEYS__Jobj for duplicates_dict in self[data_section_key][table_name]]
            else:
                temp_list = self[data_section_key][table_name]
            df = pandas.DataFrame.from_records(temp_list)
            # if table_name == 'Metabolites' and self._metabolite_header:
            #     df.columns = ['Metabolite'] + self._metabolite_header
            # elif table_name == 'Extended' and self._extended_metabolite_header:
            #     df.columns = ['Metabolite'] + self._extended_metabolite_header
            # elif table_name == 'Data':
            #     if 'BINNED' in data_section_key and self._binned_header:
            #         df.columns = ['Bin range (ppm)'] + self._binned_header
            #     elif self._samples:
            #         df.columns = ['Metabolite'] + self._samples
            return df
        return pandas.DataFrame()
    
    def get_metabolites_as_pandas(self):
        """Return the Metabolites table as a pandas.DataFrame.
        
        Note that if there are duplicate column names, they will have a string appended to the 
        end of the name like {{{_\d+_}}}.
        
        :return: The list of dicts for the Metabolites table as a pandas.DataFrame.
        :rtype: pandas.DataFrame
        """
        return self.get_table_as_pandas('Metabolites')
    
    def get_extended_as_pandas(self):
        """Return the Extended table as a pandas.DataFrame.
        
        Note that if there are duplicate column names, they will have a string appended to the 
        end of the name like {{{_\d+_}}}.
        
        :return: The list of dicts for the Extended table as a pandas.DataFrame.
        :rtype: pandas.DataFrame
        """
        return self.get_table_as_pandas('Extended')
    
    def get_metabolites_data_as_pandas(self):
        """Return the Data table as a pandas.DataFrame.
        
        Note that if there are duplicate column names, they will have a string appended to the 
        end of the name like {{{_\d+_}}}.
        
        :return: The list of dicts for the Data table as a pandas.DataFrame.
        :rtype: pandas.DataFrame
        """
        return self.get_table_as_pandas('Data')
    
    def validate(self, section_schema_mapping=section_schema_mapping, verbose=True, metabolites=True):
        """Validate the instance.
        
        :param dict section_schema_mapping: Dictionary that provides mapping between section name and schema definition.
        :param bool verbose: whether to be verbose or not.
        :param bool metabolites: whether to validate metabolites section.
        :return: Validated file and errors if verbose is False.
        :rtype: :py:class:`collections.OrderedDict`, _io.StringIO
        """
        return validate_file(
                    mwtabfile=self,
                    section_schema_mapping=section_schema_mapping,
                    verbose=verbose,
                    metabolites=metabolites
                )
    
    @classmethod
    def from_dict(cls, input_dict):
        """Create a new MWTabFile instance from input_dict.
        
        :param dict input_dict: Dictionary to create the new instance from.
        :return: New instance of MWTabFile
        :rtype: :class:`~mwtab.mwtab.MWTabFile`
        """
        new_mwtabfile = cls("Internal dictionary. ID: " + str(id(input_dict)))
        new_mwtabfile.update(input_dict)
        return new_mwtabfile
    
    def read_from_str(self, input_str):
        """Read input_str into a :class:`~mwtab.mwtab.MWTabFile` instance.

        :return: None
        :rtype: :py:obj:`None`
        """
        if not input_str:
            raise ValueError("Blank input string retrieved from source.")

        mwtab_str = self._is_mwtab(input_str)
        json_str = self._is_json(input_str)

        if not input_str:
            pass
        elif json_str:
            self.update(json_str)
        elif mwtab_str:
            self._build_mwtabfile(mwtab_str)
        else:
            raise TypeError("Unknown file format")
        
    def read(self, filehandle):
        """Read data into a :class:`~mwtab.mwtab.MWTabFile` instance.

        :param filehandle: file-like object.
        :type filehandle: :py:class:`io.TextIOWrapper`, :py:class:`gzip.GzipFile`,
                          :py:class:`bz2.BZ2File`, :py:class:`zipfile.ZipFile`
        :return: None
        :rtype: :py:obj:`None`
        """
        input_str = filehandle.read()
        self.read_from_str(input_str)

        # if not input_str:
        #     raise ValueError("Blank input string retrieved from source.")

        # mwtab_str = self._is_mwtab(input_str)
        # json_str = self._is_json(input_str)

        # if not input_str:
        #     pass
        # elif json_str:
        #     self.update(json_str)
        # elif mwtab_str:
        #     self._build_mwtabfile(mwtab_str)
        # else:
        #     raise TypeError("Unknown file format")

        filehandle.close()

    def write(self, filehandle, file_format):
        """Write :class:`~mwtab.mwtab.MWTabFile` data into file.

        :param filehandle: file-like object.
        :type filehandle: :py:class:`io.TextIOWrapper`
        :param str file_format: Format to use to write data: `mwtab` or `json`.
        :return: None
        :rtype: :py:obj:`None`
        """
        try:
            if file_format == "json":
                json_str = self._to_json()
                filehandle.write(json_str)
            elif file_format == "mwtab":
                mwtab_str = self._to_mwtab()
                filehandle.write(mwtab_str)
            else:
                raise TypeError("Unknown file format.")
        except IOError:
            raise IOError('"filehandle" parameter must be writable.')
        filehandle.close()

    def writestr(self, file_format):
        """Write :class:`~mwtab.mwtab.MWTabFile` data into string.

        :param str file_format: Format to use to write data: `mwtab` or `json`.
        :return: String representing the :class:`~mwtab.mwtab.MWTabFile` instance.
        :rtype: :py:class:`str`
        """
        if file_format == "json":
            json_str = self._to_json()
            return json_str
        elif file_format == "mwtab":
            mwtab_str = self._to_mwtab()
            return mwtab_str
        else:
            raise TypeError("Unknown file format.")

    def _build_mwtabfile(self, mwtab_str):
        """Build :class:`~mwtab.mwtab.MWTabFile` instance.

        :param mwtab_str: String in `mwtab` format.
        :type mwtab_str: :py:class:`str` or :py:class:`bytes`
        :return: instance of :class:`~mwtab.mwtab.MWTabFile`.
        :rtype: :class:`~mwtab.mwtab.MWTabFile`
        """
        mwtab_file = self
        lexer = tokenizer(mwtab_str, self._default_dict_type)
        token = next(lexer)

        while token.key != "!#ENDFILE":
            if token.key.startswith("#"):
                name = token.key[1:]
                section = self._build_block(name, lexer)
                # if section:
                if name == "METABOLITES":
                    data_section = next((n for n in mwtab_file.keys() if "METABOLITE_DATA" in n or "BINNED_DATA" in n), None)
                    if data_section:
                        for key in section.keys():
                            mwtab_file[data_section][key] = section[key]
                elif name == "NMR":
                    mwtab_file["NM"] = section
                elif name == "END":
                    pass
                else:
                    mwtab_file[name] = section
            token = next(lexer)
        return mwtab_file

    def _build_block(self, name, lexer):
        """Build individual text block of :class:`~mwtab.mwtab.MWTabFile` instance.

        :param name: name of the block, used for errors.
        :type name: str
        :param lexer: instance of the mwtab tokenizer.
        :type lexer: :func:`~mwtab.tokenizer.tokenizer`
        :return: Section dictionary.
        :rtype: :py:class:`collections.OrderedDict`
        """
        section = OrderedDict()
        token = next(lexer)
        # alias = {
        #     "subject_type": "Subject ID",
        #     "local_sample_id": "Sample ID",
        #     "factors": "Factors",
        #     "additional_sample_data": "Additional sample data",
        # }

        while not token.key.startswith("#ENDSECTION"):

            # TODO: Move to separate method (no longer works the same way as the other possibilities in _build_block())
            if token.key.startswith("SUBJECT_SAMPLE_FACTORS"):
                if type(section) != list:
                    section = list()
                # sample_dict = OrderedDict({alias[token._fields[x]]: token[x] for x in range(1, len(token._fields))})
                # if not sample_dict.get("Additional sample data"):
                #     del sample_dict["Additional sample data"]
                section.append(token.value)

            elif token.key.endswith("_START"):
                section_name = token.key
                data = []

                token = next(lexer)
                header = list(token.value)
                # Sometimes there can be extra tabs at the end of the line that results in 
                # an empty string as the token value, so remove it.
                while not header[-1]:
                    header.pop()
                metabolite_header = ["Metabolite"] + header[1:]
                if self.compatability_mode:
                    temp_header_dict = DuplicatesDict()
                    for temp_header in metabolite_header:
                        temp_header_dict[temp_header] = ''
                    metabolite_header = list(temp_header_dict.raw_keys())
                
                loop_count = 0
                while not token.key.endswith("_END"):
                    # Sometimes there can be extra tabs at the end of the line that results in 
                    # an empty string as the token value, so remove it.
                    token_value = list(token.value)
                    while not token_value[-1]:
                        token_value.pop()
                    
                    if token.key == "Bin range(ppm)" and "BINNED_DATA" in section_name and loop_count < 2:
                        self._binned_header = token_value[1:]
                        self._samples = self._binned_header
                    # Have seen Factors section in incorrect sections such as METABOLITES, 
                    # and seen multiple Factors sections in a single METABOLITE_DATA section.
                    # So just grab the one near the top.
                    elif token.key == "Factors" and "METABOLITE_DATA" in section_name and loop_count < 3:
                        self._factors = {}
                        for i, factor_string in enumerate(token_value[1:]):
                            factor_pairs = factor_string.split(" | ")
                            factor_dict = self._default_dict_type()
                            # if self.compatability_mode:
                            #     factor_dict = DuplicatesDict()
                            # else:
                            #     factor_dict = OrderedDict()
                            for pair in factor_pairs:
                                factor_key, factor_value = pair.split(":")
                                # if self.compatability_mode:
                                #     factor_dict.set(factor_key.strip(), factor_value.strip(), ordered_dict=True)
                                # else:
                                #     factor_dict[factor_key.strip()] = factor_value.strip()
                                factor_dict[factor_key.strip()] = factor_value.strip()
                            self._factors[header[i+1]] = factor_dict
                    
                    elif token.key == "Samples" and "METABOLITE_DATA" in section_name and loop_count < 3:
                        self._samples = token_value[1:]
                    
                    elif token.key.lower() == "metabolite_name" and "METABOLITES" in section_name and loop_count < 2:
                        self._metabolite_header = token_value[1:]
                    
                    elif token.key.lower() == "metabolite_name" and "EXTENDED" in section_name and loop_count < 2:
                        self._extended_metabolite_header = token_value[1:]
                        
                    else:
                        # temp_dict = self._default_dict_type()
                        # # temp_dict = DuplicatesDict() if self.compatability_mode else OrderedDict()
                        # token_len = len(token_value)
                        # for i, header_name in enumerate(metabolite_header):
                        #     # if self.compatability_mode:
                        #     #     if i < token_len:
                        #     #         temp_dict.set(header_name, token_value[i], ordered_dict=True)
                        #     #     else:
                        #     #         temp_dict.set(header_name, '', ordered_dict=True)
                        #     # else:
                        #     #     if i < token_len:
                        #     #         temp_dict[header_name] = token_value[i]
                        #     #     else:
                        #     #         temp_dict[header_name] = ""
                        #     if i < token_len:
                        #         temp_dict[header_name] = token_value[i]
                        #     else:
                        #         temp_dict[header_name] = ""
                        
                        token_len = len(token_value)
                        temp_dict = self._default_dict_type()
                        for item in zip_longest(metabolite_header, token_value, fillvalue=''):
                            temp_dict[item[0]] = item[1]
                        # temp_dict = self._default_dict_type(zip_longest(metabolite_header, token_value, fillvalue=''))
                        data.append(temp_dict)
                        
                        if token_len > len(metabolite_header):
                            self._short_headers.add(section_name)

                    token = next(lexer)
                    loop_count += 1
                
                # if self.compatability_mode:
                #     data = [DuplicatesDict(data_dict) for data_dict in data]
                if token.key.startswith("METABOLITES"):
                    section["Metabolites"] = data
                elif token.key.startswith("EXTENDED_"):
                    section["Extended"] = data
                else:
                    section["Data"] = data

            elif token.key.endswith("_RESULTS_FILE"):
                key, values = token
                results_file_dict = {}
                for i, value in enumerate(values):
                    if not value:
                        continue
                    if i == 0 and ':' not in value:
                        results_file_dict["filename"] = value
                    else:
                        split_value = value.split(':')
                        if len(split_value) == 2:
                            results_file_dict[split_value[0]] = split_value[1]
                        else:
                            results_file_dict[split_value[0]] = ":".join(split_value[1:])
                section[key] = results_file_dict

            else:
                key, value, = token
                if key in section:
                    if section[key] == value:
                        if name in self._duplicate_sub_sections:
                            self._duplicate_sub_sections[name][key] = value
                        else:
                            self._duplicate_sub_sections[name] = {key : value}
                    if name.endswith('WORKBENCH'):
                        section[key] = value
                    else:
                        section[key] += " {}".format(value)
                else:
                    section[key] = value

            # load token(s) (from parsing of next line in file)
            token = next(lexer)

        return section

    def print_file(self, f=sys.stdout, file_format="mwtab"):
        """Print :class:`~mwtab.mwtab.MWTabFile` into a file or stdout.

        :param f: writable file-like stream.
        :type f: :py:class:`io.StringIO`
        :param str file_format: Format to use: `mwtab` or `json`.
        :return: None
        :rtype: :py:obj:`None`
        """
        if file_format == "mwtab":
            for key in self:
                if key == "SUBJECT_SAMPLE_FACTORS":
                    print("#SUBJECT_SAMPLE_FACTORS:         \tSUBJECT(optional)[tab]SAMPLE[tab]FACTORS(NAME:VALUE pairs separated by |)[tab]Additional sample data", file=f)
                    self.print_subject_sample_factors(key, f=f, file_format=file_format)
                else:
                    if key == "METABOLOMICS WORKBENCH":
                        print(self.header, file=f)
                    elif key == "NM":
                        print("#NMR", file=f)
                    else:
                        print("#{}".format(key), file=f)

                    self.print_block(key, f=f, file_format=file_format)
            print("#END", file=f)

        elif file_format == "json":
            print(self._to_json(), file=f)

    def print_subject_sample_factors(self, section_key, f=sys.stdout, file_format="mwtab"):
        """Print `mwtab` `SUBJECT_SAMPLE_FACTORS` section into a file or stdout.

        :param str section_key: Section name.
        :param f: writable file-like stream.
        :type f: :py:class:`io.StringIO`
        :param str file_format: Format to use: `mwtab` or `json`.
        :return: None
        :rtype: :py:obj:`None`
        """
        if file_format == "mwtab":
            for item in self[section_key]:
                formatted_items = []
                for k in item.keys():
                    if k in ["Subject ID", "Sample ID"]:
                        formatted_items.append(str(item[k]))
                    elif k == "Factors":
                        factors = []
                        # if self.compatability_mode:
                        #     factor_dict = item[k].getObject()
                        # else:
                        #     factor_dict = item[k]
                        factor_dict = item[k]
                        for k2, value in factor_dict.items():
                            # match_re =  re.match(DUPLICATE_KEY_REGEX, k2)
                            # if match_re:
                            #     k2 = match_re.group(1)
                            
                            # if isinstance(item[k][k2], _duplicate_key_list):
                            #     for value in item[k][k2]:
                            #         factors.append("{}:{}".format(k2, value))
                            # else:
                            #     factors.append("{}:{}".format(k2, item[k][k2]))
                            factors.append("{}:{}".format(k2, value))
                        formatted_items.append(" | ".join(factors))
                    elif k == "Additional sample data":
                        additional_sample_data = []
                        # if self.compatability_mode:
                        #     add_dict = item[k].getObject()
                        # else:
                        #     add_dict = item[k]
                        add_dict = item[k]
                        for k2, value in add_dict.items():
                            # match_re =  re.match(DUPLICATE_KEY_REGEX, k2)
                            # if match_re:
                            #     k2 = match_re.group(1)
                            
                            # if isinstance(item[k][k2], _duplicate_key_list):
                            #     for value in item[k][k2]:
                            #         additional_sample_data.append("{}={}".format(k2, value))
                            # else:
                            #     additional_sample_data.append("{}={}".format(k2, item[k][k2]))
                            additional_sample_data.append("{}={}".format(k2, value))
                        formatted_items.append("; ".join(additional_sample_data))
                line = "{}{}\t{}".format(section_key, 11 * " ", "\t".join(formatted_items))
                # for file missing "Additional sample data" items
                if len(formatted_items) < 4:
                    line += "\t"
                print(line, file=f)

    def print_block(self, section_key, f=sys.stdout, file_format="mwtab"):
        """Print `mwtab` section into a file or stdout.

        :param str section_key: Section name.
        :param f: writable file-like stream.
        :type f: :py:class:`io.StringIO`
        :param str file_format: Format to use: `mwtab` or `json`.
        :return: None
        :rtype: :py:obj:`None`
        """
        if file_format == "mwtab":
            for key, value in self[section_key].items():
                if section_key == "METABOLOMICS WORKBENCH" and key not in ("VERSION", "CREATED_ON"):
                    continue

                if key in ("VERSION", "CREATED_ON"):
                    cw = 20 - len(key)
                elif key == "Units":
                    cw = 33 - len(section_key+":UNITS")
                else:
                    cw = 30 - len(key)

                if key.endswith("_RESULTS_FILE"):
                    # if isinstance(value, dict):
                    #     print("{}{}               \t{}\t{}:{}".format(self.prefixes.get(section_key, ""),
                    #                                                   *[i for pair in value.items() for i in pair]), file=f)
                    # else:
                    #     print("{}{}{}\t{}".format(self.prefixes.get(section_key, ""), key, cw * " ", value), file=f)
                    
                    results_string = self.prefixes.get(section_key, "") + key + cw * " " + "\t"
                    results_string += self._create_result_file_string(section_key, key, "mwtab")
                    
                    # split_value = value.split(':')
                    # result_value = value.split(' ')[0]
                    # for i, piece in enumerate(split_value):
                    #     for result_key in self.result_file_keys:
                    #         if piece.endswith(result_key):
                    #             if i == 0:
                    #                 result_value += "\t"
                    #             result_value += result_key + ":"
                    #             if i+2 < len(split_value):
                    #                 next_piece = split_value[i+1]
                    #                 for result_key_2 in self.result_file_keys:
                    #                     if next_piece.endswith(result_key_2):
                    #                         result_value += next_piece.replace(result_key_2, '').strip() + '\t'
                    #                         break
                    #             else:
                    #                 result_value += split_value[-1]
                    #             break
                    # results_string += result_value
                    print(results_string, file=f)

                # prints #MS_METABOLITE_DATA, #NMR_METABOLITE_DATA, or #NMR_BINNED_DATA sections
                elif key == "Units":
                    print("{}:UNITS{}\t{}".format(section_key, cw * " ", value), file=f)
                elif key == "Data":
                    print("{}_START".format(section_key), file=f)

                    if "METABOLITE" in section_key:
                        # prints "Samples" line at head of data section
                        sample_names = None
                        if self._samples is not None:
                            sample_names = self._samples
                        elif self[section_key][key]:
                            # if self.compatability_mode:
                            #     sample_names = [k for k in self[section_key][key][0].getObject().keys()][1:]
                            # else:
                            #     sample_names = [k for k in self[section_key][key][0].keys()][1:]
                            sample_names = [k for k in self[section_key][key][0].keys()][1:]
                        if sample_names:
                            print("\t".join(["Samples"] + sample_names), file=f)
                            # prints "Factors" line at head of data section
                            if self._factors is not None:
                                # factors_dict = self._factors
                                factors_dict = {}
                                for sample_name, sample_factors in self._factors.items():
                                    factors_for_sample = []
                                    # if self.compatability_mode:
                                    #     factor_dict = sample_factors.getObject()
                                    # else:
                                    #     factor_dict = sample_factors
                                    factor_dict = sample_factors
                                    for factor_key, factor_value in factor_dict.items():
                                        # match_re =  re.match(DUPLICATE_KEY_REGEX, factor_key)
                                        # if match_re:
                                        #     factor_key = match_re.group(1)
                                        factors_for_sample.append(factor_key + ":" + factor_value)
                                    factors_dict[sample_name] = ' | '.join(factors_for_sample)
                            else:
                                factors_dict = {}
                                for i in self['SUBJECT_SAMPLE_FACTORS']:
                                    sample = i['Sample ID']
                                    factors_for_sample = []
                                    # if self.compatability_mode:
                                    #     factor_dict = i['Factors'].getObject()
                                    # else:
                                    #     factor_dict = i['Factors']
                                    factor_dict = i['Factors']
                                    for factor_key, factor_value in factor_dict.items():
                                        # match_re =  re.match(DUPLICATE_KEY_REGEX, factor_key)
                                        # if match_re:
                                        #     factor_key = match_re.group(1)
                                        factors_for_sample.append(factor_key + ":" + factor_value)
                                    factors_dict[sample] = ' | '.join(factors_for_sample)
                                
                                # factors_dict = {i["Sample ID"]: i["Factors"] for i in self["SUBJECT_SAMPLE_FACTORS"]}
                            
                            factors_list = []
                            for k in sample_names:
                                # Need to make sure the factors always line up with samples.
                                if k not in factors_dict:
                                    factors_list = []
                                    break
                                # factors = [fk + ":" + factors_dict[k][fk] for fk in factors_dict[k].keys()]
                                # factors_list.append(" | ".join(factors))
                                factors_list.append(factors_dict[k])
                            if factors_list:
                                print("\t".join(["Factors"] + factors_list), file=f)
                        
                        for k, i in enumerate(self[section_key][key]):
                            # print("\t".join([i[k] for k in i.keys()]), file=f)
                            
                            # if self.compatability_mode:
                            #     print("\t".join(i.getObject().values()), file=f)
                            # else:
                            #     print("\t".join(i.values()), file=f)
                            print("\t".join(i.values()), file=f)

                    else:  # NMR_BINNED_DATA
                        # Only print if there is data to print.
                        if self._binned_header is not None:
                            binned_header = self._binned_header
                        elif self[section_key][key]:
                            # if self.compatability_mode:
                            #     binned_header = [k for k in self[section_key][key][0].getObject().keys()][1:]
                            # else:
                            #     binned_header = [k for k in self[section_key][key][0].keys()][1:]
                            binned_header = [k for k in self[section_key][key][0].keys()][1:]
                        
                        print("\t".join(["Bin range(ppm)"] + binned_header), file=f)
                        
                        for i in self[section_key][key]:
                            # print("\t".join([i[k] for k in i.keys()]), file=f)
                            
                            # if self.compatability_mode:
                            #     print("\t".join(i.getObject().values()), file=f)
                            # else:
                            #     print("\t".join(i.values()), file=f)
                            print("\t".join(i.values()), file=f)

                    print("{}_END".format(section_key), file=f)

                # prints #METABOLITES section
                elif key in ("Metabolites", "Extended"):
                    if key == "Metabolites":
                        print("#METABOLITES", file=f)
                        print("METABOLITES_START", file=f)
                    else:
                        print("EXTENDED_{}_START".format(section_key), file=f)
                    
                    if key == "Metabolites" and self._metabolite_header is not None:
                        metabolite_header = self._metabolite_header
                    elif key == "Extended" and self._extended_metabolite_header is not None:
                        metabolite_header = self._extended_metabolite_header
                    elif self[section_key][key]:
                        # if self.compatability_mode:
                        #     metabolite_header = []
                        #     for header in self[section_key][key][0].getObject().keys():
                        #         match_re =  re.match(DUPLICATE_KEY_REGEX, header)
                        #         if match_re:
                        #             header = match_re.group(1)
                        #         metabolite_header.append(header)
                        #     metabolite_header = metabolite_header[1:]
                        # else:
                        #     metabolite_header = [k for k in self[section_key][key][0].keys()][1:]
                        metabolite_header = [k for k in self[section_key][key][0].keys()][1:]
                    else:
                        metabolite_header = []
                    print("\t".join(["metabolite_name"] + metabolite_header), file=f)
                    
                    for i in self[section_key][key]:
                        # print("\t".join(i[k] for k in i.keys()), file=f)
                        
                        # if self.compatability_mode:
                        #     print("\t".join(i.getObject().values()), file=f)
                        # else:
                        #     print("\t".join(i.values()), file=f)
                        print("\t".join(i.values()), file=f)

                    if key == "Metabolites":
                        print("METABOLITES_END", file=f)
                    else:
                        print("EXTENDED_{}_END".format(section_key), file=f)

                else:
                    # Filenames don't get split.
                    if len(str(value)) > 80 and not key.endswith("_FILENAME"):
                        words = str(value).split(" ")
                        length = 0
                        line = list()
                        for word in words:
                            if length + len(word) + len(line) - 1 < 80:
                                line.append(word)
                                length += len(word)
                            else:
                                # Long filenames were printing 2 lines, with the top one being blank.
                                # I fixed this by adding a check to skip filenames, but just in case I also don't 
                                # let empty lines be printed.
                                if line:
                                    print("{}{}{}\t{}".format(self.prefixes.get(section_key, ""), key, cw * " ", " ".join(line)), file=f)
                                line = [word]
                                length = len(word)
                        print("{}{}{}\t{}".format(self.prefixes.get(section_key, ""), key, cw * " ", " ".join(line)),
                              file=f)
                    else:
                        print("{}{}{}\t{}".format(self.prefixes.get(section_key, ""), key, cw * " ", value), file=f)

        elif file_format == "json":
            print(json.dumps(self[section_key], sort_keys=SORT_KEYS, indent=INDENT), file=f)

    def _to_json(self):
        """Save :class:`~mwtab.mwtab.MWTabFile` into JSON string.

        :return: JSON string.
        :rtype: :py:class:`str`
        """
        self._set_key_order()
        temp = copy.deepcopy(OrderedDict(self))
        # Result files ends up being a dictionary, but needs to printed as a string.
        for section_key, section_value in temp.items():
            if isinstance(section_value, (dict, OrderedDict)):
                for key in section_value:
                    if key.endswith("_RESULTS_FILE"):
                        temp[section_key][key] = self._create_result_file_string(section_key, key, "json")
        
        if self.compatability_mode:
            # Handle duplicate keys in SSF.
            temp = jdks.JSON_DUPLICATE_KEYS(temp)
            # for ssf_dict in temp.getObject()["SUBJECT_SAMPLE_FACTORS"]:
                
            #     if "Factors" in ssf_dict:
            #         factor_dict = ssf_dict["Factors"]
            #         new_factor_dict = DuplicatesDict()
            #         for key, value in factor_dict.items():
            #             if isinstance(value, _duplicate_key_list):
            #                 for item in value:
            #                     new_factor_dict.set(key, item, ordered_dict=True)
            #             else:
            #                 new_factor_dict.set(key, value, ordered_dict=True)
            #             ssf_dict["Factors"] = new_factor_dict
                
            #     if "Additional sample data" in ssf_dict:
            #         additional_dict = ssf_dict["Additional sample data"]
            #         new_additional_dict = DuplicatesDict()
            #         for key, value in additional_dict.items():
            #             if isinstance(value, _duplicate_key_list):
            #                 for item in value:
            #                     new_additional_dict.set(key, item, ordered_dict=True)
            #             else:
            #                 new_additional_dict.set(key, value, ordered_dict=True)
            #             ssf_dict["Additional sample data"] = new_additional_dict
            
            json_string = temp.dumps(sort_keys=SORT_KEYS, indent=INDENT, default=_JSON_serializer_for_dupe_class)
            json_string = re.sub(r'("Additional sample data"): "\{(.*)\}"', _match_process, json_string)
            json_string = re.sub(r'("Factors"): "\{(.*)\}"', _match_process, json_string)
            return json_string
        else:
            return json.dumps(temp, sort_keys=SORT_KEYS, indent=INDENT)
    

    def _to_mwtab(self):
        """Save :class:`~mwtab.mwtab.MWTabFile` in `mwtab` formatted string.

        :return: NMR-STAR string.
        :rtype: :py:class:`str`
        """
        self._set_key_order()
        mwtab_str = io.StringIO()
        self.print_file(mwtab_str)
        return mwtab_str.getvalue()

    @staticmethod
    def _is_mwtab(string):
        """Test if input string is in `mwtab` format.

        :param string: Input string.
        :type string: :py:class:`str` or :py:class:`bytes`
        :return: Input string if in mwTab format or False otherwise.
        :rtype: :py:class:`str` or :py:obj:`False`
        """
        if isinstance(string, str):
            lines = string.replace("\r", "\n").split("\n")
        elif isinstance(string, bytes):
            lines = string.decode("utf-8").replace("\r", "\n").split("\n")
        else:
            raise TypeError("Expecting <class 'str'> or <class 'bytes'>, but {} was passed".format(type(string)))

        lines = [line for line in lines if line]
        header = lines[0]

        if header.startswith("#METABOLOMICS WORKBENCH"):
            return "\n".join(lines)
        return False

    @staticmethod
    def _is_json(string, compatability_mode=False):
        """Test if input string is in JSON format.

        :param string: Input string.
        :type string: :py:class:`str` or :py:class:`bytes`
        :param compatability_mode: if true, replace some dictionaries with DuplicatesDict.
        :type text: py:class:`bool`
        :return: Input string if in JSON format or False otherwise.
        :rtype: :py:class:`str` or :py:obj:`False`
        """
        try:
            if isinstance(string, bytes):
                if compatability_mode:
                    json_str = json.loads(string.decode("utf-8"), object_pairs_hook=_handle_duplicate_keys)
                else:
                    json_str = json.loads(string.decode("utf-8"))
            elif isinstance(string, str):
                if compatability_mode:
                    json_str = json.loads(string, object_pairs_hook=_handle_duplicate_keys)
                else:
                    json_str = json.loads(string)
            else:
                raise TypeError("Expecting <class 'str'> or <class 'bytes'>, but {} was passed".format(type(string)))
            
            if compatability_mode:
                for i, ssf_dict in enumerate(json_str['SUBJECT_SAMPLE_FACTORS']):
                    if not isinstance(ssf_dict['Factors'], DuplicatesDict):
                        json_str['SUBJECT_SAMPLE_FACTORS'][i]['Factors'] = DuplicatesDict(ssf_dict['Factors'])
                    if 'Additional sample data' in ssf_dict and not isinstance(ssf_dict['Additional sample data'], DuplicatesDict):
                        json_str['SUBJECT_SAMPLE_FACTORS'][i]['Additional sample data'] = DuplicatesDict(ssf_dict['Additional sample data'])
                
                data_key = [n for n in json_str.keys() if "METABOLITE_DATA" in n or "BINNED_DATA" in n]
                if data_key:
                    data_key = data_key[0]
                    for section_key, section_value in json_str[data_key].items():
                        if section_key in ["Data", "Metabolites", "Extended"]:
                            for i, data_dict in enumerate(section_value):
                                json_str[data_key][section_key][i] = DuplicatesDict(json_str[data_key][section_key][i])
                
            
            return json_str
        except ValueError:
            return False
    
    def _create_result_file_string(self, section_key, key, to_format):
        """
        """
        delimiter = "\t" if to_format == "mwtab" else " "
        
        if "filename" in self[section_key][key]:
            new_value = self[section_key][key]["filename"]
        else:
            new_value = ""
        
        other_pair_strings = []
        for result_key, result_value in self[section_key][key].items():
            if result_key != "filename":
                other_pair_strings.append(result_key + ":" + result_value)
        
        joined_pairs = delimiter.join(other_pair_strings)
        
        if new_value:
            if joined_pairs:
                new_value += delimiter + joined_pairs
        elif joined_pairs:
            new_value += joined_pairs
        return new_value
    
    def _set_key_order(self):
        """
        """
        # This is a way to build the key_order diction from the schema, but it can't be used 
        # until we only support Python > 3.6 because we have to be able to rely on ordered 
        # dict behavior.
        # from .mwschema import section_schema_mapping
        # key_order = {}
        # for key, schema in section_schema_mapping.items():
        #     json_schema = schema.json_schema("asdf")
        #     if "properties" in json_schema:
        #         key_order[key] = {}
        #         for json_property, property_dict in json_schema["properties"].items():
        #             if "items" in property_dict:
        #                 key_order[key][json_property] = list(property_dict["items"]["properties"].keys())
        #             else:
        #                 key_order[key][json_property] = []
        #     elif "items" in json_schema:
        #         key_order[key] = {json_property:[] for json_property in json_schema["items"]["properties"]}
               
        key_order = \
            {'METABOLOMICS WORKBENCH': {},
             'PROJECT': {},
             'STUDY': {},
             'SUBJECT': {},
             'SUBJECT_SAMPLE_FACTORS': {'Subject ID': [],
              'Sample ID': [],
              'Factors': [],
              'Additional sample data': []},
             'COLLECTION': {},
             'TREATMENT': {},
             'SAMPLEPREP': {},
             'CHROMATOGRAPHY': {},
             'ANALYSIS': {},
             'MS': {},
             'NM': {},
             'MS_METABOLITE_DATA': {'Units': [],
              'Data': ['Metabolite', 'Bin range(ppm)'],
              'Metabolites': ['Metabolite', 'Bin range(ppm)'],
              'Extended': ['Metabolite']},
             'NMR_METABOLITE_DATA': {'Units': [],
              'Data': ['Metabolite', 'Bin range(ppm)'],
              'Metabolites': ['Metabolite', 'Bin range(ppm)'],
              'Extended': ['Metabolite']},
             'NMR_BINNED_DATA': {'Units': [], 'Data': ['Metabolite', 'Bin range(ppm)']}}
    
        for key, sub_keys in key_order.items():
            if key in self:
                # SUBJECT_SAMPLE_FACTORS is the only list as of now.
                if isinstance(self[key], list):
                    temp_list = []
                    for i, element in enumerate(self[key]):
                        temp_list.append(OrderedDict())
                        for sub_key in sub_keys:
                            if sub_key in element:
                                temp_list[i][sub_key] = element[sub_key]
                    del self[key]
                    self[key] = temp_list
                else:
                    temp_dict = OrderedDict()
                    for sub_key, sub_sub_keys in sub_keys.items():
                        if sub_key in self[key]:
                            # For the Data, Metabolites, and Extended sections.
                            if isinstance(self[key][sub_key], list):
                                temp_list = []
                                for i, element in enumerate(self[key][sub_key]):
                                    temp_list.append(OrderedDict())
                                    for sub_sub_key in sub_sub_keys:
                                        # if self.compatability_mode:
                                        #     element_dict = element.getObject()
                                        # else:
                                        #     element_dict = element
                                        if sub_sub_key in element:
                                            temp_list[i][sub_sub_key] = element[sub_sub_key]
                                    # Add the elements that aren't in the key_order.
                                    for unordered_key in element:
                                        if unordered_key not in temp_list[i]:
                                            temp_list[i][unordered_key] = element[unordered_key]
                                    if self.compatability_mode:
                                        temp_list[i] = DuplicatesDict(temp_list[i])
                                temp_dict[sub_key] = temp_list
                            else:
                                temp_dict[sub_key] = self[key][sub_key]
                    # Add any unknown sub_keys that aren't in key_order.
                    for unordered_sub_key in self[key]:
                        if unordered_sub_key not in temp_dict:
                            temp_dict[unordered_sub_key] = self[key][unordered_sub_key]
                            
                    del self[key]
                    self[key] = temp_dict
        
        # Handle unrecognized keys, put them at the end.
        keys_to_move = []
        for key in self:
            if key not in key_order:
                keys_to_move.append(key)
        
        for key in keys_to_move:
            temp = self[key]
            del self[key]
            self[key] = temp
    
    
    
    
    
