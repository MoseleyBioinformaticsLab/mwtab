#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
mwtab.mwtab
~~~~~~~~~~~

This module provides the :class:`~mwtab.mwtab.MWTabFile` class
that stores the data from a single ``mwTab`` formatted file in the
form of an :py:class:`dict`. Data can be accessed
directly from the :class:`~mwtab.mwtab.MWTabFile` instance using
bracket accessors.

The data is divided into a series of "sections" which each contain a
number of "key-value"-like pairs. Also, the file contains a specially
formatted ``SUBJECT_SAMPLE_FACTOR`` block and blocks of data between 
``*_START`` and ``*_END``.
"""

from __future__ import print_function, division, unicode_literals
import io
import sys
import json
import re
import copy
from itertools import zip_longest

import pandas

from .tokenizer import tokenizer, _results_file_line_to_dict
from .validator import validate_file
from .mwschema import ms_required_schema, nmr_required_schema
from .duplicates_dict import DuplicatesDict, DUPLICATE_KEY_REGEX

SORT_KEYS = False
INDENT = 4


# From https://stackoverflow.com/questions/14902299/json-loads-allows-duplicate-keys-in-a-dictionary-overwriting-the-first-value
def _handle_duplicate_keys(ordered_pairs):
    """Object pairs hook for the json package to be able to read in duplicate keys for dictionaries."""
    d = DuplicatesDict()
    s = set()
    for k, v in ordered_pairs:
        d[k] = v
        s.add(k)
    
    if len(s) == len(ordered_pairs):
        return d.data
    else:
        return d


def _parse_header_input(input_str):
    """Attempt to parse a header string into a dict."""
    match_re = re.match(r"#METABOLOMICS WORKBENCH( )?([^: ]+ )?([A-Z_]+:\w+ ?)*", input_str)
    if match_re:
        key_values = re.findall(r" (\w*:\w*)", input_str)
        temp_dict = {}
        for key_value in key_values:
            key, new_value = key_value.split(":")
            temp_dict[key] = new_value
        
        return temp_dict
    else:
        raise ValueError(r"Header cannot be set because it is not of the form \"#METABOLOMICS WORKBENCH( )?([^: ]+ )?([A-Z_]+:\w+ ?)*\"")



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
                if isinstance(obj["METABOLOMICS WORKBENCH"], dict):
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
            


class MWTabFile(dict):
    """MWTabFile class that stores data from a single ``mwTab`` formatted file in
    the form of a dictionary.
    
    Parameters:
        source: A string that should be the file path to the mwtab file that will be read in.
        duplicate_keys: If True, use a special dictionary type that can handle duplicate keys. 
                        If you are uisng this class to build an mwtab file by hand, don't set 
                        this to True. This was added because some files already upload to the 
                        Metabolomics Workbench erroneously have duplicate keys and the class 
                        needed to be able to read write them back out correctly.
        force: If True, replace non-dictionary values in METABOLITES_DATA, METABOLITES, and EXTENDED 
               tables with empty dicts on JSON read in.
    
    Attributes:
        source: A string that should be the file path to the mwtab file that was read in.
        study_id: A managed property. The study ID is stored in the METABOLOMICS WORKBENCH 
                  key of this class and the JSON version of an mwTab file. This property 
                  is provided as a convenience to access the study ID.
        analysis_id: A managed property. The analysis ID is stored in the METABOLOMICS WORKBENCH 
                     key of this class and the JSON version of an mwTab file. This property 
                     is provided as a convenience to access the analysis ID.
        header: A managed property. It is provided as a convenience to be able to view the 
                header line of the file. You can also set the METABOLOMICS WORKBENCH key by 
                setting this property. It will parse the string you assign into the dictionary 
                that belongs in the METABOLOMICS WORKBENCH key. Assuming the provided string is 
                a correctly generated header line.
        data_section_key: A simple property that will give you the key to the data section. 
                          Either one of "MS_METABOLITE_DATA", "NMR_METABOLITE_DATA", or "NMR_BINNED_DATA", 
                          or None if none of those were found.
    
    Special Notes:
        In general this class has the same structure as mwTab JSON, but there are a few exceptions to that.
        One is that the _RESULTS_FILE subsection, whether in the MS or NM section, is a dictionary in this 
        class with keys for the elements found. In the mwTab JSON this is just a string. The possible 
        keys the dictionary could have are: "filename", "UNITS", "Has m/z", "Has RT", and "RT units". 
        If the _RESULTS_FILE line did not have these keys, then they won't be in the dictionary.
        
        The Metabolomics Workbench has deprecated mwTab files with NMR_BINNED_DATA sections, but for 
        the few that do exist if you read them in using this class, the dictionaries in the ['NMR_BINNED_DATA']['Data'] 
        list of dicts will have keys for both "Metabolite" and "Bin range(ppm)". This is because the 
        JSON version from the Metabolomics Workbench uses "Bin range(ppm)" and not "Metabolite", but 
        we wanted to present a seemless unified interface for this class regardless of the analysis type. 
        They print out with only the "Bin range(ppm)" keys to match what the Metabolomics Workbench provides, 
        but internally both keys will be there. If they somehow become different, you will see a message 
        about it when you try to write the file out.
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
    table_names = ['Data', 'Extended', 'Metabolites']
    
    study_id = MWTabProperty()
    analysis_id = MWTabProperty()
    header = MWTabProperty()

    def __init__(self, source, duplicate_keys=False, force=False, *args, **kwds):
        """File initializer.

        :param str source: Source a `MWTabFile` instance was created from.
        """
        super(MWTabFile, self).__init__(*args, **kwds)
        self.source = source
        self._force = force
        self._factors = None
        self._samples = None
        self._raw_samples = None
        self._metabolite_header = None
        self._raw_metabolite_header = None
        self._extended_metabolite_header = None
        self._raw_extended_metabolite_header = None
        self._binned_header = None
        self._raw_binned_header = None
        self._short_headers = set()
        self._duplicate_sub_sections = {}
        self._input_format = 'json'
        self._duplicate_keys = duplicate_keys
        if duplicate_keys:
            self._default_dict_type = DuplicatesDict
        else:
            self._default_dict_type = dict
    
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
    
    def set_table_from_pandas(self, df, table_name, clear_header=False):
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
        self.set_table_from_pandas(df, 'Metabolites', clear_header)
    
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
        self.set_table_from_pandas(df, 'Extended', clear_header)
    
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
        self.set_table_from_pandas(df, 'Data', clear_header)
    
    def get_table_as_pandas(self, table_name):
        """Return the given table_name as a pandas.DataFrame.
        
        table_name must be one of "Metabolites", "Extended", or "Data". Note that 
        if there are duplicate column names, they will have a string appended to the 
        end of the name like {{{_\\d+_}}}.
        
        :param str table_name: the name of the table to return as a pandas.DataFrame.
        :return: The list of dicts for the given table_name as a pandas.DataFrame.
        :rtype: pandas.DataFrame
        """
        data_section_key = self.data_section_key
        if data_section_key and table_name in self[data_section_key]:
            if self._duplicate_keys:
                temp_list = [duplicates_dict.data for duplicates_dict in self[data_section_key][table_name]]
            else:
                temp_list = self[data_section_key][table_name]
            df = pandas.DataFrame.from_records(temp_list)
            return df
        return pandas.DataFrame()
    
    def get_metabolites_as_pandas(self):
        """Return the Metabolites table as a pandas.DataFrame.
        
        Note that if there are duplicate column names, they will have a string appended to the 
        end of the name like {{{_\\d+_}}}.
        
        :return: The list of dicts for the Metabolites table as a pandas.DataFrame.
        :rtype: pandas.DataFrame
        """
        return self.get_table_as_pandas('Metabolites')
    
    def get_extended_as_pandas(self):
        """Return the Extended table as a pandas.DataFrame.
        
        Note that if there are duplicate column names, they will have a string appended to the 
        end of the name like {{{_\\d+_}}}.
        
        :return: The list of dicts for the Extended table as a pandas.DataFrame.
        :rtype: pandas.DataFrame
        """
        return self.get_table_as_pandas('Extended')
    
    def get_metabolites_data_as_pandas(self):
        """Return the Data table as a pandas.DataFrame.
        
        Note that if there are duplicate column names, they will have a string appended to the 
        end of the name like {{{_\\d+_}}}.
        
        :return: The list of dicts for the Data table as a pandas.DataFrame.
        :rtype: pandas.DataFrame
        """
        return self.get_table_as_pandas('Data')
    
    def validate(self, ms_schema: dict = ms_required_schema, nmr_schema: dict = nmr_required_schema, verbose: bool = True) -> (str, list[dict]):
        """Validate the instance.
        
        Args:
            ms_schema: jsonschema to validate both the base parts of the file and the MS specific parts of the file.
            nmr_schema: jsonschema to validate both the base parts of the file and the NMR specific parts of the file.
            verbose: whether to be verbose or not.
        
        Returns: 
            Error messages as a single string and error messages in JSON form. If verbose is True, then the single string will be None.
        """
        return validate_file(
                    mwtabfile=self,
                    ms_schema = ms_schema,
                    nmr_schema = nmr_schema,
                    verbose = verbose,
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
        self._input_format = 'mwtab' if mwtab_str else 'json'
        json_str = self._is_json(input_str, self._duplicate_keys, self._force)

        if json_str:
            self.update(json_str)
            
            metabolite_header = [column if not column.endswith('}}}') else re.match(DUPLICATE_KEY_REGEX, column).group(1) 
                                 for column in self.get_metabolites_as_pandas().columns]
            self._raw_metabolite_header = metabolite_header if metabolite_header else None
            self._metabolite_header = metabolite_header[1:] if metabolite_header[1:] else None
            
            extended_metabolite_header = [column if not column.endswith('}}}') else re.match(DUPLICATE_KEY_REGEX, column).group(1) 
                                          for column in self.get_extended_as_pandas().columns]
            self._raw_extended_metabolite_header = extended_metabolite_header if extended_metabolite_header else None
            self._extended_metabolite_header = extended_metabolite_header[1:] if extended_metabolite_header[1:] else None
            
            samples = [column if not column.endswith('}}}') else re.match(DUPLICATE_KEY_REGEX, column).group(1) 
                       for column in self.get_metabolites_data_as_pandas().columns]
            self._raw_samples = samples if samples else None
            self._samples = samples[1:] if samples[1:] else None
            
            if (data_section_key := self.data_section_key) and "BINNED" in data_section_key:
                self._binned_header = self._samples
                self._raw_binned_header = self._raw_samples
                
                for i in range(len(self['NMR_BINNED_DATA']['Data'])):
                    if 'Bin range(ppm)' in self['NMR_BINNED_DATA']['Data'][i]:
                        self['NMR_BINNED_DATA']['Data'][i]['Metabolite'] = self['NMR_BINNED_DATA']['Data'][i]['Bin range(ppm)']
        
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
        
        # Sometimes the results file line is in the wrong spot, look for it in DATA and move it if so.
        # Looked at changing the tokenizer to use the two letter code on the line, 'MS:RESULTS_FILE', 
        # so we could put it in the correct spot when building, but there are other examples of two 
        # letter codes that are in the right spot, but have the wrong 2 letter code (AN000012).
        if self.data_section_key:
            results_file_key = None
            for key in self[self.data_section_key]:
                if key.endswith('_RESULTS_FILE'):
                    results_file_key = key
                    if 'MS' in self:
                        section_key = 'MS'
                    elif 'NM' in self:
                        section_key = 'NM'
                    
                    temp = self[self.data_section_key][key]
                    self[section_key][key] = temp
            if results_file_key:
                del self[self.data_section_key][results_file_key]
                
        return mwtab_file

    def _build_block(self, name, lexer):
        """Build individual text block of :class:`~mwtab.mwtab.MWTabFile` instance.

        :param name: name of the block, used for errors.
        :type name: str
        :param lexer: instance of the mwtab tokenizer.
        :type lexer: :func:`~mwtab.tokenizer.tokenizer`
        :return: Section dictionary.
        :rtype: :py:class:`dict`
        """
        section = {}
        token = next(lexer)
        
        if 'SUBJECT_SAMPLE_FACTORS' in self:
            ssf_samples = [value['Sample ID'] for value in self['SUBJECT_SAMPLE_FACTORS']]
        
        while not token.key.startswith("#ENDSECTION"):

            if token.key.startswith("SUBJECT_SAMPLE_FACTORS"):
                if type(section) != list:
                    section = list()
                section.append(token.value)

            elif token.key.endswith("_START"):
                section_name = token.key[0:-6]
                data = []

                token = next(lexer)
                header = list(token.value)
                # Sometimes there can be extra tabs at the end of the line that results in 
                # an empty string as the token value, so remove it.
                while not header[-1]:
                    header.pop()
                metabolite_header = ["Metabolite"] + header[1:]
                
                loop_count = 0
                while not token.key.endswith("_END"):
                    # Sometimes there can be extra tabs at the end of the line that results in 
                    # an empty string as the token value, so remove it.
                    token_value = list(token.value)
                    while not token_value[-1]:
                        token_value.pop()
                    
                    is_header = False
                    if "BINNED_DATA" in section_name and 'EXTENDED' not in section_name and loop_count < 2:
                        if loop_count < 1:
                            self._raw_binned_headers = token_value
                            self._raw_samples = self._raw_binned_headers
                        if token.key == "Bin range(ppm)":
                            # self._binned_header = token_value[1:]
                            # self._samples = self._binned_header
                            is_header = True
                    # Have seen Factors section in incorrect sections such as METABOLITES, 
                    # and seen multiple Factors sections in a single METABOLITE_DATA section.
                    # So just grab the one near the top.
                    elif token.key == "Factors" and "METABOLITE_DATA" in section_name and loop_count < 2:
                        self._factors = {}
                        for i, factor_string in enumerate(token_value[1:]):
                            factor_pairs = factor_string.split("| ")
                            factor_dict = self._default_dict_type()
                            for pair in factor_pairs:
                                factor_key, factor_value = pair.split(":")
                                factor_dict[factor_key.strip()] = factor_value.strip()
                            self._factors[header[i+1]] = factor_dict
                        is_header = True
                    
                    elif "METABOLITE_DATA" in section_name and 'EXTENDED' not in section_name and loop_count < 2:
                        if loop_count < 1:
                            self._raw_samples = token_value
                        # The last check for len(token_value) == 1 is for ones like AN000788.
                        if (any(sample in ssf_samples for sample in token_value[1:]) or \
                           (len(token_value) == 1 and token_value[0] in ['Samples', 'metabolite name', 'metabolite_name']) or \
                           (len(ssf_samples) == 0 and token_value[0] in ['Samples', 'metabolite name', 'metabolite_name'])):
                            # self._samples = token_value[1:]
                            is_header = True
                    
                    elif "METABOLITES" in section_name and loop_count < 2:
                        if loop_count < 1:
                            self._raw_metabolite_header = token_value
                        if token.key.lower() == "metabolite_name":
                            # self._metabolite_header = token_value[1:]
                            is_header = True
                    
                    elif "EXTENDED" in section_name and loop_count < 2:
                        if loop_count < 1:
                            self._raw_extended_metabolite_header = token_value
                        if token.key.lower() == "metabolite_name":
                            # self._extended_metabolite_header = token_value[1:]
                            is_header = True
                    
                    
                    
                    # TODO remove.
                    # if token.key == "Bin range(ppm)" and "BINNED_DATA" in section_name and loop_count < 2:
                    #     self._binned_header = token_value[1:]
                    #     self._samples = self._binned_header
                    # # Have seen Factors section in incorrect sections such as METABOLITES, 
                    # # and seen multiple Factors sections in a single METABOLITE_DATA section.
                    # # So just grab the one near the top.
                    # elif token.key == "Factors" and "METABOLITE_DATA" in section_name and loop_count < 2:
                    #     self._factors = {}
                    #     for i, factor_string in enumerate(token_value[1:]):
                    #         factor_pairs = factor_string.split("| ")
                    #         factor_dict = self._default_dict_type()
                    #         for pair in factor_pairs:
                    #             factor_key, factor_value = pair.split(":")
                    #             factor_dict[factor_key.strip()] = factor_value.strip()
                    #         self._factors[header[i+1]] = factor_dict
                    
                    # # The last check for len(token_value) == 1 is for ones like AN000788.
                    # elif "METABOLITE_DATA" in section_name and 'EXTENDED' not in section_name and self._samples is None and \
                    #      loop_count < 2 and (any(sample in ssf_samples for sample in token_value[1:]) or \
                    #      (len(token_value) == 1 and token_value[0] in ['Samples', 'metabolite name', 'metabolite_name'])):
                    #     self._samples = token_value[1:]
                    
                    # elif token.key.lower() == "metabolite_name" and "METABOLITES" in section_name and loop_count < 2:
                    #     self._metabolite_header = token_value[1:]
                    
                    # elif token.key.lower() == "metabolite_name" and "EXTENDED" in section_name and loop_count < 2:
                    #     self._extended_metabolite_header = token_value[1:]
                        
                    if not is_header:                        
                        token_len = len(token_value)
                        temp_dict = self._default_dict_type()
                        for item in zip_longest(metabolite_header, token_value, fillvalue=''):
                            temp_dict[item[0]] = item[1]
                        data.append(temp_dict)
                        
                        if token_len > len(metabolite_header):
                            self._short_headers.add(section_name)

                    token = next(lexer)
                    loop_count += 1
                
                # This makes it so all dicitonaries have the same number of values.
                # Let's say row 3 looks like {'Metabolite': 'asdf', 'col1': 'qwer', '': 2345}
                # The rows above don't have the '' entry, this code makes it so they do.
                if self._duplicate_keys:
                    data = [duplicates_dict.data for duplicates_dict in data]
                data_df = pandas.DataFrame.from_records(data).fillna('').astype(str)
                data = data_df.to_dict(orient='records')
                if self._duplicate_keys:
                    data = [DuplicatesDict(data_dict) for data_dict in data]
                min_header = [column if not column.endswith('}}}') else re.match(DUPLICATE_KEY_REGEX, column).group(1) 
                              for column in data_df.columns]
                min_header = min_header[1:] if min_header[1:] else None
                
                if token.key.startswith("METABOLITES"):
                    section["Metabolites"] = data
                    self._metabolite_header = min_header
                elif token.key.startswith("EXTENDED_"):
                    section["Extended"] = data
                    self._extended_metabolite_header = min_header
                else:
                    section["Data"] = data
                    self._samples = min_header
                    if "BINNED_DATA" in section_name:
                        self._binned_header = min_header

            elif token.key.endswith("_RESULTS_FILE"):
                key, results_file_dict = token
                section[key] = results_file_dict

            else:
                key, value = token
                if key in section:
                    if section[key] == value:
                        self._duplicate_sub_sections.setdefault(name, {})
                        self._duplicate_sub_sections[name][key] = value
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
                    
                    if isinstance(self[key], dict):
                        self.print_block(key, f=f, file_format=file_format)
                    else:
                        raise TypeError(f'Key/section "{key}" is not a dictionary. It cannot be translated to the mwTab format.')
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
                        factor_dict = item[k]
                        for k2, value in factor_dict.items():
                            factors.append("{}:{}".format(k2, value))
                        formatted_items.append(" | ".join(factors))
                    elif k == "Additional sample data":
                        additional_sample_data = []
                        add_dict = item[k]
                        for k2, value in add_dict.items():
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
                    
                    results_string = self.prefixes.get(section_key, "") + key + cw * " " + "\t"
                    results_string += self._create_result_file_string(section_key, key, "mwtab")
                    print(results_string, file=f)

                # prints #MS_METABOLITE_DATA, #NMR_METABOLITE_DATA, or #NMR_BINNED_DATA sections
                elif key == "Units":
                    print("{}:UNITS{}\t{}".format(section_key, cw * " ", value), file=f)
                elif key == "Data":
                    print("{}_START".format(section_key), file=f)

                    if "METABOLITE" in section_key:
                        # prints "Samples" line at head of data section
                        sample_names = []
                        if self._samples is not None:
                            sample_names = self._samples
                        elif self[section_key][key]:
                            sample_names = [k for k in self[section_key][key][0].keys()][1:]
                        print("\t".join(["Samples"] + sample_names), file=f)
                        if sample_names:
                            # prints "Factors" line at head of data section
                            if self._factors is not None:
                                factors_dict = {}
                                for sample_name, sample_factors in self._factors.items():
                                    factors_for_sample = []
                                    factor_dict = sample_factors
                                    for factor_key, factor_value in factor_dict.items():
                                        factors_for_sample.append(factor_key + ":" + factor_value)
                                    factors_dict[sample_name] = ' | '.join(factors_for_sample)
                            else:
                                factors_dict = {}
                                for i in self['SUBJECT_SAMPLE_FACTORS']:
                                    sample = i['Sample ID']
                                    factors_for_sample = []
                                    factor_dict = i['Factors']
                                    for factor_key, factor_value in factor_dict.items():
                                        factors_for_sample.append(factor_key + ":" + factor_value)
                                    factors_dict[sample] = ' | '.join(factors_for_sample)
                            
                            factors_list = []
                            for k in sample_names:
                                # Need to make sure the factors always line up with samples.
                                if k not in factors_dict:
                                    factors_list = []
                                    break
                                factors_list.append(factors_dict[k])
                            if factors_list:
                                print("\t".join(["Factors"] + factors_list), file=f)
                        
                        for k, i in enumerate(self[section_key][key]):
                            print("\t".join(i.values()), file=f)

                    else:  # NMR_BINNED_DATA
                        # Only print if there is data to print.
                        if self._binned_header is not None:
                            binned_header = self._binned_header
                        elif self[section_key][key]:
                            binned_header = [k for k in self[section_key][key][0].keys()][1:]
                        
                        print("\t".join(["Bin range(ppm)"] + binned_header), file=f)
                        
                        for i in self[section_key][key]:
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
                        metabolite_header = [k for k in self[section_key][key][0].keys()][1:]
                    else:
                        metabolite_header = []
                    print("\t".join(["metabolite_name"] + metabolite_header), file=f)
                    
                    for i in self[section_key][key]:
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
        
        # Note that indent cannot be None or json will use a version of the 
        # encoder written in C and DuplicatesDict will not be printed correctly.
        elif file_format == "json":
            print(json.dumps(self[section_key], sort_keys=SORT_KEYS, indent=INDENT), file=f)

    def _to_json(self):
        """Save :class:`~mwtab.mwtab.MWTabFile` into JSON string.

        :return: JSON string.
        :rtype: :py:class:`str`
        """
        self._set_key_order()
        temp = copy.deepcopy(self)
        # Result files ends up being a dictionary, but needs to printed as a string.
        for section_key, section_value in temp.items():
            if isinstance(section_value, dict):
                for key in section_value:
                    if key.endswith("_RESULTS_FILE"):
                        temp[section_key][key] = self._create_result_file_string(section_key, key, "json")
        
        if 'NMR_BINNED_DATA' in temp:
            for i in range(len(temp['NMR_BINNED_DATA']['Data'])):
                if 'Bin range(ppm)' in temp['NMR_BINNED_DATA']['Data'][i]:
                    if temp['NMR_BINNED_DATA']['Data'][i]['Metabolite'] != temp['NMR_BINNED_DATA']['Data'][i]['Bin range(ppm)']:
                        print("Warning: The \"Metabolite\" key and \"Bin range(ppm)\" "
                              f"key in ['NMR_BINNED_DATA']['Data'][{i}] are different "
                              "values. Only the value in \"Bin range(ppm)\" will be written out.")
                    del temp['NMR_BINNED_DATA']['Data'][i]['Metabolite']
                else:
                    new_dict = self._default_dict_type()
                    new_dict['Bin range(ppm)'] = temp['NMR_BINNED_DATA']['Data'][i]['Metabolite']
                    del temp['NMR_BINNED_DATA']['Data'][i]['Metabolite']
                    # Note that the update method cannot be used here because it can mess up DuplicatesDict.
                    for key, value in temp['NMR_BINNED_DATA']['Data'][i].items():
                        new_dict[key] = value
                    temp['NMR_BINNED_DATA']['Data'][i] = new_dict
        
        # Note that indent cannot be None or json will use a version of the 
        # encoder written in C and DuplicatesDict will not be printed correctly.
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
    def _is_json(string, duplicate_keys=False, force=False):
        """Test if input string is in JSON format.

        :param string: Input string.
        :type string: :py:class:`str` or :py:class:`bytes`
        :param duplicate_keys: if true, replace some dictionaries with DuplicatesDict.
        :type text: py:class:`bool`
        :param force: if true, replace non-dictionary values in METABOLITES_DATA, METABOLITES, and EXTENDED with empty dicts.
        :type text: py:class:`bool`
        :return: Input string if in JSON format or False otherwise.
        :rtype: :py:class:`str` or :py:obj:`False`
        """
        try:
            if isinstance(string, bytes):
                if duplicate_keys:
                    json_str = json.loads(string.decode("utf-8"), object_pairs_hook=_handle_duplicate_keys)
                else:
                    json_str = json.loads(string.decode("utf-8"))
            elif isinstance(string, str):
                if duplicate_keys:
                    json_str = json.loads(string, object_pairs_hook=_handle_duplicate_keys)
                else:
                    json_str = json.loads(string)
            else:
                raise TypeError("Expecting <class 'str'> or <class 'bytes'>, but {} was passed".format(type(string)))
            
            if duplicate_keys and isinstance(json_str, DuplicatesDict):
                raise KeyError('Given JSON contains duplicate keys at the highest level.')
            
            if duplicate_keys:
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
                        indexes_to_delete = []
                        for i, data_dict in enumerate(section_value):
                            if isinstance(data_dict, dict):
                                if duplicate_keys and not isinstance(data_dict, DuplicatesDict):
                                    json_str[data_key][section_key][i] = DuplicatesDict(json_str[data_key][section_key][i])
                            else:
                                if force:
                                    indexes_to_delete.append(i)
                                else:
                                    raise TypeError(f'Given JSON contains non-dictionary values in ["{data_key}"]["{section_key}"].')
                        for index in indexes_to_delete:
                            del json_str[data_key][section_key][index]
            
            # Have to tokenize the results file.
            # "ST000071_AN000111_Results.txt UNITS:Peak area Has m/z:Yes Has RT:Yes RT units:Minutes"
            for section_key, section_value in json_str.items():
                if isinstance(section_value, dict):
                    results_file_keys = [key for key in section_value if key.endswith('_RESULTS_FILE')]
                    if results_file_keys:
                        results_file_key = results_file_keys[0]
                        results_file_value = json_str[section_key][results_file_key]
                        results_file_dict = _results_file_line_to_dict(results_file_value)
                        
                        json_str[section_key][results_file_key] = results_file_dict
                
            return json_str
        except ValueError:
            return False
    
    def _create_result_file_string(self, section_key: str, key: str, to_format: str) -> str:
        """Build result file string from its dictionary parts.
        
        Args:
            section_key: The section where the result file dicitonary is.
            key: The key where the result file is.
            to_format: If 'mwtab', then formats the string for mwTab, otherwise formats for JSON.
        
        Returns:
            A string to go in either the mwTab or JSON file.
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
        """Set the key order to a specific order.
        
        Sets the key order to a certain order for better reproducibility and consistency.
        """               
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
                        temp_list.append({})
                        for sub_key in sub_keys:
                            if sub_key in element:
                                temp_list[i][sub_key] = element[sub_key]
                    del self[key]
                    self[key] = temp_list
                else:
                    temp_dict = {}
                    for sub_key, sub_sub_keys in sub_keys.items():
                        if sub_key in self[key]:
                            # For the Data, Metabolites, and Extended sections.
                            if isinstance(self[key][sub_key], list):
                                temp_list = []
                                for i, element in enumerate(self[key][sub_key]):
                                    temp_list.append({})
                                    for sub_sub_key in sub_sub_keys:
                                        if sub_sub_key in element:
                                            temp_list[i][sub_sub_key] = element[sub_sub_key]
                                    # Add the elements that aren't in the key_order.
                                    for unordered_key in element:
                                        if unordered_key not in temp_list[i]:
                                            temp_list[i][unordered_key] = element[unordered_key]
                                    if self._duplicate_keys:
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
    
    def __deepcopy__(self, memo):
        new_tabfile = MWTabFile(self.source, self._duplicate_keys)
        memo[id(new_tabfile)] = new_tabfile
        for key, value in self.items():
            new_tabfile[key] = copy.deepcopy(value, memo)
        for key, value in self.__dict__.items():
            new_tabfile.__dict__[key] = copy.deepcopy(value, memo)
        return new_tabfile
    
    def __copy__(self):
        new_tabfile = MWTabFile(self.source, self._duplicate_keys)
        for key, value in self.items():
            new_tabfile[key] = copy.copy(value)
        for key, value in self.__dict__.items():
            new_tabfile.__dict__[key] = copy.copy(value)
        return new_tabfile
    
    
    
    
