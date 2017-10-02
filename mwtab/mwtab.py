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

from .tokenizer import tokenizer


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
        "MS_METABOLITE_DATA": "",
        "NMR_BINNED_DATA": "",
        "METABOLITES": ""
    }

    def __init__(self, source, *args, **kwds):
        """File initializer.

        :param str source: Source a `MWTabFile` instance was created from.
        """
        super(MWTabFile, self).__init__(*args, **kwds)
        self.source = source
        self.study_id = ""
        self.analysis_id = ""
        self.header = ""

    def read(self, filehandle):
        """Read data into a :class:`~mwtab.mwtab.MWTabFile` instance.

        :param filehandle: file-like object.
        :type filehandle: :py:class:`io.TextIOWrapper`, :py:class:`gzip.GzipFile`,
                          :py:class:`bz2.BZ2File`, :py:class:`zipfile.ZipFile`
        :return: None
        :rtype: :py:obj:`None`
        """
        input_str = filehandle.read()
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

        self.study_id = self["METABOLOMICS WORKBENCH"].get("STUDY_ID")
        self.analysis_id = self["METABOLOMICS WORKBENCH"].get("ANALYSIS_ID")
        self.header = self["METABOLOMICS WORKBENCH"].get("HEADER")

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
        lexer = tokenizer(mwtab_str)
        token = next(lexer)

        while token.key != "!#ENDFILE":
            if token.key.startswith("#"):
                name = token.key[1:]
                section = self._build_block(lexer)
                if section:
                    mwtab_file[name] = section
            token = next(lexer)
        return mwtab_file

    def _build_block(self, lexer):
        """Build individual text block of :class:`~mwtab.mwtab.MWTabFile` instance.

        :param lexer: instance of the mwtab tokenizer.
        :type lexer: :func:`~mwtab.tokenizer.tokenizer`
        :return: Section dictionary.
        :rtype: :py:class:`collections.OrderedDict`
        """
        section = OrderedDict()
        token = next(lexer)

        while not token.key.startswith("#ENDSECTION"):

            if token.key.startswith("SUBJECT_SAMPLE_FACTORS"):
                header = token._fields[1:]
                values = token[1:]
                section.setdefault("SUBJECT_SAMPLE_FACTORS", [])
                section["SUBJECT_SAMPLE_FACTORS"].append(OrderedDict(zip(header, values)))

            elif token.key.endswith("_START"):
                section_key = token.key
                section.setdefault(section_key, OrderedDict())
                data = []

                token = next(lexer)
                header = list(token.value)

                while not token.key.endswith("_END"):
                    if token.key == "Samples":
                        header[0] = "metabolite_name"
                        section[section_key][token.key] = list(token.value[1:])
                    elif token.key == "Factors":
                        section[section_key][token.key] = list(token.value[1:])
                    elif token.key in ("metabolite_name", "Bin range(ppm)"):
                        section[section_key]["Fields"] = list(token.value)
                    else:
                        data.append(OrderedDict(zip(header, token.value)))

                    token = next(lexer)
                section[section_key]["DATA"] = data

            elif token.key.endswith("_RESULTS_FILE"):
                if len(token) == 4:
                    key, value, extrakey, extravalue = token
                    section[key] = OrderedDict([(key, value), (extrakey,extravalue)])
                else:
                    key, value = token
                    section[key] = value

            else:
                key, value, = token
                if key in section:
                    section[key] += "\n{}".format(value)
                else:
                    section[key] = value
            token = next(lexer)
        return section

    def print_file(self, f=sys.stdout, file_format="mwtab"):
        """Print :class:`~mwtab.mwtab.MWTabFile` into a file or stdout.

        :param io.StringIO f: writable file-like stream.
        :param str file_format: Format to use: `mwtab` or `json`.
        :param f: Print to file or stdout.
        :param int tw: Tab width.
        :return: None
        :rtype: :py:obj:`None`
        """
        if file_format == "mwtab":
            for key in self:
                if key == "SUBJECT_SAMPLE_FACTORS":
                    print("#SUBJECT_SAMPLE_FACTORS:         \tSUBJECT(optional)[tab]SAMPLE[tab]FACTORS(NAME:VALUE pairs separated by |)[tab]Additional sample data", file=f)
                elif key == "METABOLOMICS WORKBENCH":
                    print(self.header, file=f)
                else:
                    print("#{}".format(key), file=f)

                self.print_block(key, f=f, file_format=file_format)
            print("#END", file=f)

        elif file_format == "json":
            print(self._to_json(), file=f)

    def print_block(self, section_key, f=sys.stdout, file_format="mwtab"):
        """Print `mwtab` section into a file or stdout.

        :param str section_key: Section name.
        :param io.StringIO f: writable file-like stream.
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
                elif key in ("SUBJECT_SAMPLE_FACTORS", ):
                    cw = 33 - len(key)
                else:
                    cw = 30 - len(key)

                if "\n" in value:
                    for line in value.split("\n"):
                        print("{}{}{}\t{}".format(self.prefixes.get(section_key, ""), key, cw * " ", line), file=f)

                elif key == "SUBJECT_SAMPLE_FACTORS":
                    for factor in value:
                        print("{}{}\t{}".format(key, cw * " ", "\t".join(factor.values())), file=f)

                elif key.endswith(":UNITS"):
                    print("{}\t{}".format(key, value), file=f)

                elif key.endswith("_RESULTS_FILE"):
                    if isinstance(value, dict):
                        print("{}{}               \t{}\t{}:{}".format(self.prefixes.get(section_key, ""),
                                                                      *[i for pair in value.items() for i in pair]), file=f)
                    else:
                        print("{}{}{}\t{}".format(self.prefixes.get(section_key, ""), key, cw * " ", value), file=f)

                elif key.endswith("_START"):
                    start_key = key
                    end_key = "{}{}".format(start_key[:-5], "END")
                    print(start_key, file=f)

                    for data_key in value:
                        if data_key in ("Samples", "Factors"):
                            print("{}\t{}".format(data_key, "\t".join(self[section_key][key][data_key])), file=f)

                        elif data_key in ("Fields", ):
                            print("{}".format("\t".join(self[section_key][key][data_key])), file=f)

                        elif data_key == "DATA":
                            for data in self[section_key][key][data_key]:
                                print("\t".join(data.values()), file=f)

                    print(end_key, file=f)
                else:
                    print("{}{}{}\t{}".format(self.prefixes.get(section_key, ""), key, cw * " ", value), file=f)

        elif file_format == "json":
            print(json.dumps(self[section_key], sort_keys=False, indent=4), file=f)

    def _to_json(self):
        """Save :class:`~mwtab.mwtab.MWTabFile` into JSON string.

        :return: JSON string.
        :rtype: :py:class:`str`
        """
        return json.dumps(self, sort_keys=False, indent=4)

    def _to_mwtab(self):
        """Save :class:`~mwtab.mwtab.MWTabFile` in `mwtab` formatted string.

        :return: NMR-STAR string.
        :rtype: :py:class:`str`
        """
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
            lines = string.split("\n")
        elif isinstance(string, bytes):
            lines = string.decode("utf-8").split("\n")
        else:
            raise TypeError("Expecting <class 'str'> or <class 'bytes'>, but {} was passed".format(type(string)))

        lines = [line for line in lines if line]
        header = lines[0]

        if header.startswith("#METABOLOMICS WORKBENCH"):
            return "\n".join(lines)
        return False

    @staticmethod
    def _is_json(string):
        """Test if input string is in JSON format.

        :param string: Input string.
        :type string: :py:class:`str` or :py:class:`bytes`
        :return: Input string if in JSON format or False otherwise.
        :rtype: :py:class:`str` or :py:obj:`False`
        """
        try:
            if isinstance(string, bytes):
                json_str = json.loads(string.decode("utf-8"), object_pairs_hook=OrderedDict)
            elif isinstance(string, str):
                json_str = json.loads(string, object_pairs_hook=OrderedDict)
            else:
                raise TypeError("Expecting <class 'str'> or <class 'bytes'>, but {} was passed".format(type(string)))
            return json_str
        except ValueError:
            return False
