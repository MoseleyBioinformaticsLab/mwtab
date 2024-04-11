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
        self._set_key_order()
        return json.dumps(self, sort_keys=False, indent=4)

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
            {'METABOLOMICS WORKBENCH': {'VERSION': [],
              'CREATED_ON': [],
              'STUDY_ID': [],
              'ANALYSIS_ID': [],
              'PROJECT_ID': [],
              'HEADER': [],
              'DATATRACK_ID': []},
             'PROJECT': {'PROJECT_TITLE': [],
              'PROJECT_TYPE': [],
              'PROJECT_SUMMARY': [],
              'INSTITUTE': [],
              'DEPARTMENT': [],
              'LABORATORY': [],
              'LAST_NAME': [],
              'FIRST_NAME': [],
              'ADDRESS': [],
              'EMAIL': [],
              'PHONE': [],
              'FUNDING_SOURCE': [],
              'PROJECT_COMMENTS': [],
              'PUBLICATIONS': [],
              'CONTRIBUTORS': [],
              'DOI': []},
             'STUDY': {'STUDY_TITLE': [],
              'STUDY_TYPE': [],
              'STUDY_SUMMARY': [],
              'INSTITUTE': [],
              'DEPARTMENT': [],
              'LABORATORY': [],
              'LAST_NAME': [],
              'FIRST_NAME': [],
              'ADDRESS': [],
              'EMAIL': [],
              'PHONE': [],
              'NUM_GROUPS': [],
              'TOTAL_SUBJECTS': [],
              'NUM_MALES': [],
              'NUM_FEMALES': [],
              'STUDY_COMMENTS': [],
              'PUBLICATIONS': [],
              'SUBMIT_DATE': []},
             'ANALYSIS': {'ANALYSIS_TYPE': [],
              'LABORATORY_NAME': [],
              'OPERATOR_NAME': [],
              'DETECTOR_TYPE': [],
              'SOFTWARE_VERSION': [],
              'ACQUISITION_DATE': [],
              'ANALYSIS_PROTOCOL_FILE': [],
              'ACQUISITION_PARAMETERS_FILE': [],
              'PROCESSING_PARAMETERS_FILE': [],
              'DATA_FORMAT': [],
              'ACQUISITION_ID': [],
              'ACQUISITION_TIME': [],
              'ANALYSIS_COMMENTS': [],
              'ANALYSIS_DISPLAY': [],
              'INSTRUMENT_NAME': [],
              'INSTRUMENT_PARAMETERS_FILE': [],
              'NUM_FACTORS': [],
              'NUM_METABOLITES': [],
              'PROCESSED_FILE': [],
              'RANDOMIZATION_ORDER': [],
              'RAW_FILE': []},
             'SUBJECT': {'SUBJECT_TYPE': [],
              'SUBJECT_SPECIES': [],
              'TAXONOMY_ID': [],
              'GENOTYPE_STRAIN': [],
              'AGE_OR_AGE_RANGE': [],
              'WEIGHT_OR_WEIGHT_RANGE': [],
              'HEIGHT_OR_HEIGHT_RANGE': [],
              'GENDER': [],
              'HUMAN_RACE': [],
              'HUMAN_ETHNICITY': [],
              'HUMAN_TRIAL_TYPE': [],
              'HUMAN_LIFESTYLE_FACTORS': [],
              'HUMAN_MEDICATIONS': [],
              'HUMAN_PRESCRIPTION_OTC': [],
              'HUMAN_SMOKING_STATUS': [],
              'HUMAN_ALCOHOL_DRUG_USE': [],
              'HUMAN_NUTRITION': [],
              'HUMAN_INCLUSION_CRITERIA': [],
              'HUMAN_EXCLUSION_CRITERIA': [],
              'ANIMAL_ANIMAL_SUPPLIER': [],
              'ANIMAL_HOUSING': [],
              'ANIMAL_LIGHT_CYCLE': [],
              'ANIMAL_FEED': [],
              'ANIMAL_WATER': [],
              'ANIMAL_INCLUSION_CRITERIA': [],
              'CELL_BIOSOURCE_OR_SUPPLIER': [],
              'CELL_STRAIN_DETAILS': [],
              'SUBJECT_COMMENTS': [],
              'CELL_PRIMARY_IMMORTALIZED': [],
              'CELL_PASSAGE_NUMBER': [],
              'CELL_COUNTS': [],
              'SPECIES_GROUP': []},
             'SUBJECT_SAMPLE_FACTORS': {'Subject ID': [],
              'Sample ID': [],
              'Factors': [],
              'Additional sample data': []},
             'COLLECTION': {'COLLECTION_SUMMARY': [],
              'COLLECTION_PROTOCOL_ID': [],
              'COLLECTION_PROTOCOL_FILENAME': [],
              'COLLECTION_PROTOCOL_COMMENTS': [],
              'SAMPLE_TYPE': [],
              'COLLECTION_METHOD': [],
              'COLLECTION_LOCATION': [],
              'COLLECTION_FREQUENCY': [],
              'COLLECTION_DURATION': [],
              'COLLECTION_TIME': [],
              'VOLUMEORAMOUNT_COLLECTED': [],
              'STORAGE_CONDITIONS': [],
              'COLLECTION_VIALS': [],
              'STORAGE_VIALS': [],
              'COLLECTION_TUBE_TEMP': [],
              'ADDITIVES': [],
              'BLOOD_SERUM_OR_PLASMA': [],
              'TISSUE_CELL_IDENTIFICATION': [],
              'TISSUE_CELL_QUANTITY_TAKEN': []},
             'TREATMENT': {'TREATMENT_SUMMARY': [],
              'TREATMENT_PROTOCOL_ID': [],
              'TREATMENT_PROTOCOL_FILENAME': [],
              'TREATMENT_PROTOCOL_COMMENTS': [],
              'TREATMENT': [],
              'TREATMENT_COMPOUND': [],
              'TREATMENT_ROUTE': [],
              'TREATMENT_DOSE': [],
              'TREATMENT_DOSEVOLUME': [],
              'TREATMENT_DOSEDURATION': [],
              'TREATMENT_VEHICLE': [],
              'ANIMAL_VET_TREATMENTS': [],
              'ANIMAL_ANESTHESIA': [],
              'ANIMAL_ACCLIMATION_DURATION': [],
              'ANIMAL_FASTING': [],
              'ANIMAL_ENDP_EUTHANASIA': [],
              'ANIMAL_ENDP_TISSUE_COLL_LIST': [],
              'ANIMAL_ENDP_TISSUE_PROC_METHOD': [],
              'ANIMAL_ENDP_CLINICAL_SIGNS': [],
              'HUMAN_FASTING': [],
              'HUMAN_ENDP_CLINICAL_SIGNS': [],
              'CELL_STORAGE': [],
              'CELL_GROWTH_CONTAINER': [],
              'CELL_GROWTH_CONFIG': [],
              'CELL_GROWTH_RATE': [],
              'CELL_INOC_PROC': [],
              'CELL_MEDIA': [],
              'CELL_ENVIR_COND': [],
              'CELL_HARVESTING': [],
              'PLANT_GROWTH_SUPPORT': [],
              'PLANT_GROWTH_LOCATION': [],
              'PLANT_PLOT_DESIGN': [],
              'PLANT_LIGHT_PERIOD': [],
              'PLANT_HUMIDITY': [],
              'PLANT_TEMP': [],
              'PLANT_WATERING_REGIME': [],
              'PLANT_NUTRITIONAL_REGIME': [],
              'PLANT_ESTAB_DATE': [],
              'PLANT_HARVEST_DATE': [],
              'PLANT_GROWTH_STAGE': [],
              'PLANT_METAB_QUENCH_METHOD': [],
              'PLANT_HARVEST_METHOD': [],
              'PLANT_STORAGE': [],
              'CELL_PCT_CONFLUENCE': [],
              'CELL_MEDIA_LASTCHANGED': []},
             'SAMPLEPREP': {'SAMPLEPREP_SUMMARY': [],
              'SAMPLEPREP_PROTOCOL_ID': [],
              'SAMPLEPREP_PROTOCOL_FILENAME': [],
              'SAMPLEPREP_PROTOCOL_COMMENTS': [],
              'PROCESSING_METHOD': [],
              'PROCESSING_STORAGE_CONDITIONS': [],
              'EXTRACTION_METHOD': [],
              'EXTRACT_CONCENTRATION_DILUTION': [],
              'EXTRACT_ENRICHMENT': [],
              'EXTRACT_CLEANUP': [],
              'EXTRACT_STORAGE': [],
              'SAMPLE_RESUSPENSION': [],
              'SAMPLE_DERIVATIZATION': [],
              'SAMPLE_SPIKING': [],
              'ORGAN': [],
              'ORGAN_SPECIFICATION': [],
              'CELL_TYPE': [],
              'SUBCELLULAR_LOCATION': []},
             'CHROMATOGRAPHY': {'CHROMATOGRAPHY_SUMMARY': [],
              'CHROMATOGRAPHY_TYPE': [],
              'INSTRUMENT_NAME': [],
              'COLUMN_NAME': [],
              'FLOW_GRADIENT': [],
              'FLOW_RATE': [],
              'COLUMN_TEMPERATURE': [],
              'METHODS_FILENAME': [],
              'SOLVENT_A': [],
              'SOLVENT_B': [],
              'METHODS_ID': [],
              'COLUMN_PRESSURE': [],
              'INJECTION_TEMPERATURE': [],
              'INTERNAL_STANDARD': [],
              'INTERNAL_STANDARD_MT': [],
              'RETENTION_INDEX': [],
              'RETENTION_TIME': [],
              'SAMPLE_INJECTION': [],
              'SAMPLING_CONE': [],
              'ANALYTICAL_TIME': [],
              'CAPILLARY_VOLTAGE': [],
              'MIGRATION_TIME': [],
              'OVEN_TEMPERATURE': [],
              'PRECONDITIONING': [],
              'RUNNING_BUFFER': [],
              'RUNNING_VOLTAGE': [],
              'SHEATH_LIQUID': [],
              'TIME_PROGRAM': [],
              'TRANSFERLINE_TEMPERATURE': [],
              'WASHING_BUFFER': [],
              'WEAK_WASH_SOLVENT_NAME': [],
              'WEAK_WASH_VOLUME': [],
              'STRONG_WASH_SOLVENT_NAME': [],
              'STRONG_WASH_VOLUME': [],
              'TARGET_SAMPLE_TEMPERATURE': [],
              'SAMPLE_LOOP_SIZE': [],
              'SAMPLE_SYRINGE_SIZE': [],
              'RANDOMIZATION_ORDER': [],
              'CHROMATOGRAPHY_COMMENTS': []},
             'MS': {'INSTRUMENT_NAME': [],
              'INSTRUMENT_TYPE': [],
              'MS_TYPE': [],
              'ION_MODE': [],
              'MS_COMMENTS': [],
              'CAPILLARY_TEMPERATURE': [],
              'CAPILLARY_VOLTAGE': [],
              'COLLISION_ENERGY': [],
              'COLLISION_GAS': [],
              'DRY_GAS_FLOW': [],
              'DRY_GAS_TEMP': [],
              'FRAGMENT_VOLTAGE': [],
              'FRAGMENTATION_METHOD': [],
              'GAS_PRESSURE': [],
              'HELIUM_FLOW': [],
              'ION_SOURCE_TEMPERATURE': [],
              'ION_SPRAY_VOLTAGE': [],
              'IONIZATION': [],
              'IONIZATION_ENERGY': [],
              'IONIZATION_POTENTIAL': [],
              'MASS_ACCURACY': [],
              'PRECURSOR_TYPE': [],
              'REAGENT_GAS': [],
              'SOURCE_TEMPERATURE': [],
              'SPRAY_VOLTAGE': [],
              'ACTIVATION_PARAMETER': [],
              'ACTIVATION_TIME': [],
              'ATOM_GUN_CURRENT': [],
              'AUTOMATIC_GAIN_CONTROL': [],
              'BOMBARDMENT': [],
              'CDL_SIDE_OCTOPOLES_BIAS_VOLTAGE': [],
              'CDL_TEMPERATURE': [],
              'DATAFORMAT': [],
              'DESOLVATION_GAS_FLOW': [],
              'DESOLVATION_TEMPERATURE': [],
              'INTERFACE_VOLTAGE': [],
              'IT_SIDE_OCTOPOLES_BIAS_VOLTAGE': [],
              'LASER': [],
              'MATRIX': [],
              'NEBULIZER': [],
              'OCTPOLE_VOLTAGE': [],
              'PROBE_TIP': [],
              'RESOLUTION_SETTING': [],
              'SAMPLE_DRIPPING': [],
              'SCAN_RANGE_MOVERZ': [],
              'SCANNING': [],
              'SCANNING_CYCLE': [],
              'SCANNING_RANGE': [],
              'SKIMMER_VOLTAGE': [],
              'TUBE_LENS_VOLTAGE': [],
              'MS_RESULTS_FILE': []},
             'NM': {'INSTRUMENT_NAME': [],
              'INSTRUMENT_TYPE': [],
              'NMR_EXPERIMENT_TYPE': [],
              'NMR_COMMENTS': [],
              'FIELD_FREQUENCY_LOCK': [],
              'STANDARD_CONCENTRATION': [],
              'SPECTROMETER_FREQUENCY': [],
              'NMR_PROBE': [],
              'NMR_SOLVENT': [],
              'NMR_TUBE_SIZE': [],
              'SHIMMING_METHOD': [],
              'PULSE_SEQUENCE': [],
              'WATER_SUPPRESSION': [],
              'PULSE_WIDTH': [],
              'POWER_LEVEL': [],
              'RECEIVER_GAIN': [],
              'OFFSET_FREQUENCY': [],
              'PRESATURATION_POWER_LEVEL': [],
              'CHEMICAL_SHIFT_REF_CPD': [],
              'TEMPERATURE': [],
              'NUMBER_OF_SCANS': [],
              'DUMMY_SCANS': [],
              'ACQUISITION_TIME': [],
              'RELAXATION_DELAY': [],
              'SPECTRAL_WIDTH': [],
              'NUM_DATA_POINTS_ACQUIRED': [],
              'REAL_DATA_POINTS': [],
              'LINE_BROADENING': [],
              'ZERO_FILLING': [],
              'APODIZATION': [],
              'BASELINE_CORRECTION_METHOD': [],
              'CHEMICAL_SHIFT_REF_STD': [],
              'BINNED_INCREMENT': [],
              'BINNED_DATA_NORMALIZATION_METHOD': [],
              'BINNED_DATA_PROTOCOL_FILE': [],
              'BINNED_DATA_CHEMICAL_SHIFT_RANGE': [],
              'BINNED_DATA_EXCLUDED_RANGE': [],
              'NMR_RESULTS_FILE': []},
             'MS_METABOLITE_DATA': {'Units': [],
              'Data': ['Metabolite', 'Bin range(ppm)'],
              'Metabolites': ['Metabolite', 'Bin range(ppm)'],
              'Extended': ['Metabolite', 'sample_id']},
             'NMR_METABOLITE_DATA': {'Units': [],
              'Data': ['Metabolite', 'Bin range(ppm)'],
              'Metabolites': ['Metabolite', 'Bin range(ppm)'],
              'Extended': ['Metabolite', 'sample_id']},
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
                                        if sub_sub_key in element:
                                            temp_list[i][sub_sub_key] = element[sub_sub_key]
                                    # Add the elements that aren't in the key_order.
                                    for unordered_key in element:
                                        if unordered_key not in temp_list[i]:
                                            temp_list[i][unordered_key] = element[unordered_key]
                                temp_dict[sub_key] = temp_list
                            else:
                                temp_dict[sub_key] = self[key][sub_key]
                    # Add any unknown sub_keys that aren't in key_order.
                    for unordered_sub_key in self[key]:
                        if unordered_sub_key not in temp_dict:
                            temp_dict[unordered_sub_key] = self[key][unordered_sub_key]
                            
                    del self[key]
                    self[key] = temp_dict
    
    
    
    
    
