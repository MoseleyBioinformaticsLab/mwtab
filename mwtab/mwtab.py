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

        try:
            self.study_id = self["METABOLOMICS WORKBENCH"].get("STUDY_ID")
            self.analysis_id = self["METABOLOMICS WORKBENCH"].get("ANALYSIS_ID")
            # self.header = self["METABOLOMICS WORKBENCH"].get("HEADER")
            self.header = " ".join(
                ["#METABOLOMICS WORKBENCH"]
                + [item[0] + ":" + item[1] for item in self["METABOLOMICS WORKBENCH"].items() if item[0] not in ["VERSION", "CREATED_ON"]]
            )
        except KeyError as e:
            raise KeyError("File missing header information \"METABOLOMICS WORKBENCH\"", e)

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
                    if name == "METABOLITES":
                        data_section = next((n for n in mwtab_file.keys() if n in ("MS_METABOLITE_DATA", "NMR_METABOLITE_DATA", "NMR_BINNED_DATA")), None)
                        if data_section:
                            for key in section.keys():
                                mwtab_file[data_section][key] = section[key]
                    elif name == "NMR":
                        mwtab_file["NM"] = section
                    else:
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
        alias = {
            "subject_type": "Subject ID",
            "local_sample_id": "Sample ID",
            "factors": "Factors",
            "additional_sample_data": "Additional sample data",
        }

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
                data = []

                token = next(lexer)
                header = list(token.value)

                while not token.key.endswith("_END"):
                    if token.key in ("Samples", "Factors", "metabolite_name", "Bin range(ppm)"):
                        pass
                    else:
                        data.append(OrderedDict(zip(["Metabolite"] + header[1:], token.value)))

                    token = next(lexer)

                if token.key.startswith("METABOLITES"):
                    section["Metabolites"] = data
                elif token.key.startswith("EXTENDED_"):
                    section["Extended"] = data
                else:
                    section["Data"] = data

            elif token.key.endswith("_RESULTS_FILE"):
                if len(token) > 2:
                    key, value, extra = token
                    section[key] = OrderedDict([(key, value)])
                    for pair in extra:
                        section[key].update((pair,))
                else:
                    key, value = token
                    section[key] = value

            else:
                key, value, = token
                if key in section:
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
                        for k2 in item[k]:
                            factors.append("{}:{}".format(k2, item[k][k2]))
                        formatted_items.append(" | ".join(factors))
                    elif k == "Additional sample data":
                        additional_sample_data = []
                        for k2 in item[k]:
                            additional_sample_data.append("{}={}".format(k2, item[k][k2]))
                        formatted_items.append(";".join(additional_sample_data))
                line = "{}{}\t{}".format(section_key, 33 * " ", "\t".join(formatted_items))
                # for file missing "Additional sample data" items
                if len(formatted_items) < 4:
                    line += "\t"
                else:
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
                    if isinstance(value, dict):
                        print("{}{}               \t{}\t{}:{}".format(self.prefixes.get(section_key, ""),
                                                                      *[i for pair in value.items() for i in pair]), file=f)
                    else:
                        print("{}{}{}\t{}".format(self.prefixes.get(section_key, ""), key, cw * " ", value), file=f)

                # prints #MS_METABOLITE_DATA, #NMR_METABOLITE_DATA, or #NMR_BINNED_DATA sections
                elif key == "Units":
                    print("{}:UNITS{}\t{}".format(section_key, cw * " ", value), file=f)
                elif key == "Data":
                    print("{}_START".format(section_key), file=f)

                    if "METABOLITE" in section_key:
                        # prints "Samples" line at head of data section
                        print("\t".join(["Samples"] + [k for k in self[section_key][key][0].keys()][1:]), file=f)
                        # prints "Factors" line at head of data section
                        factors_list = ["Factors"]
                        factors_dict = {i["Sample ID"]: i["Factors"] for i in self["SUBJECT_SAMPLE_FACTORS"]}
                        for k in [k for k in self[section_key][key][0].keys()][1:]:
                            factors = [fk + ":" + factors_dict[k][fk] for fk in factors_dict[k].keys()]
                            factors_list.append(" | ".join(factors))
                        print("\t".join(factors_list), file=f)
                        for i in self[section_key][key]:
                            print("\t".join([i[k] for k in i.keys()]), file=f)

                    else:  # NMR_BINNED_DATA
                        print("\t".join(["Bin range(ppm)"] + [k for k in self[section_key][key][0].keys()][1:]), file=f)
                        for i in self[section_key][key]:
                            print("\t".join([i[k] for k in i.keys()]), file=f)

                    print("{}_END".format(section_key), file=f)

                # prints #METABOLITES section
                elif key in ("Metabolites", "Extended"):
                    if key == "Metabolites":
                        print("#METABOLITES", file=f)
                        print("METABOLITES_START", file=f)
                    else:
                        print("EXTENDED_{}_START".format(section_key), file=f)

                    print("\t".join(["metabolite_name"] + [k for k in self[section_key][key][0].keys()][1:]), file=f)
                    for i in self[section_key][key]:
                        print("\t".join(i[k] for k in i.keys()), file=f)

                    if key == "Metabolites":
                        print("METABOLITES_END", file=f)
                    else:
                        print("EXTENDED_{}_END".format(section_key), file=f)

                else:
                    if len(str(value)) > 80:
                        words = str(value).split(" ")
                        length = 0
                        line = list()
                        for word in words:
                            if length + len(word) + len(line) - 1 <= 80:
                                line.append(word)
                                length += len(word)
                            else:
                                print("{}{}{}\t{}".format(self.prefixes.get(section_key, ""), key, cw * " ", " ".join(line)), file=f)
                                line = [word]
                                length = len(word)
                        print("{}{}{}\t{}".format(self.prefixes.get(section_key, ""), key, cw * " ", " ".join(line)),
                              file=f)
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
