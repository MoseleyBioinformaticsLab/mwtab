#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# TODO: Add qualifications for output items per specific context input item.
"""
mwtab.mwrest
~~~~~~~~~~~~

This module provides routines for accessing the Metabolomics Workbench REST API.

See https://www.metabolomicsworkbench.org/tools/MWRestAPIv1.0.pdf for details.
"""

from collections import OrderedDict
from . import fileio
import re
import json


VERBOSE = False
BASE_URL = "https://www.metabolomicsworkbench.org/rest/"


def analysis_ids(base_url=BASE_URL):
    """
    Method for retrieving a list of analysis ids for every current analysis in Metabolomics Workbench.

    :param str base_url: Base url to Metabolomics Workbench REST API.
    :return: List of every available Metabolomics Workbench analysis identifier.
    :rtype: :py:class:`list`
    """
    st_an_dict = _pull_study_analysis(base_url)
    analyses = list()
    [analyses.extend(st_an_dict[k]) for k in st_an_dict.keys()]

    if VERBOSE:
        print("Found {} Analysis Files to be Downloaded".format(len(analyses)))

    return analyses


def study_ids(base_url=BASE_URL):
    """
    Method for retrieving a list of study ids for every current study in Metabolomics Workbench.

    :param str base_url: Base url to Metabolomics Workbench REST API.
    :return: List of every available Metabolomics Workbench study identifier.
    :rtype: :py:class:`list`
    """
    st_an_dict = _pull_study_analysis(base_url)
    studies = list(st_an_dict.keys())

    if VERBOSE:
        print("Found {} Study Files to be Downloaded".format(len(studies)))

    return studies


def _pull_study_analysis(base_url=BASE_URL):
    """
    Method for requesting a JSON string containing all study ids and analysis ids from Metabolomics Workbench's REST
    API. Requests a JSON file which contains a list of studies and their accompanying analyses. The JSON file is
    converted into a python object (dict) which can then be parsed to create a dictionary with the form study id (key):
    analysis id(s) (values).

    :param str base_url: Base url to Metabolomics Workbench REST API.
    :return: Dictionary of study ids (keys) and lists of analyses (value).
    :rtype: :py:class:`dict`
    """
    url = GenericMWURL(
        {"context": "study", "input_item": "study_id", "input_value": "ST", "output_item": "analysis"},
        base_url
    ).url
    mwrestfile = next(fileio.read_mwrest(url, **{"convertJSON": True}))
    json_object = mwrestfile._is_json(mwrestfile.text)

    study_analysis_dict = dict()
    for k in json_object.keys():
        if study_analysis_dict.get(json_object[k]["study_id"]):
            study_analysis_dict[json_object[k]["study_id"]].append(json_object[k]["analysis_id"])
        else:
            study_analysis_dict[json_object[k]["study_id"]] = [json_object[k]["analysis_id"]]

    return study_analysis_dict


def generate_mwtab_urls(input_items, base_url=BASE_URL, output_format='txt'):
    """
    Method for generating URLS to be used to retrieve `mwtab` files for analyses and
    studies through the REST API of the Metabolomics Workbench database.

    :param list input_items: List of Metabolomics Workbench input values for mwTab files.
    :param str base_url: Base url to Metabolomics Workbench REST API.
    :param str output_format: Output format for the mwTab files to be retrieved in.
    :return: Metabolomics Workbench REST URL string(s).
    :rtype: :py:class:`str`
    """
    for input_item in input_items:
        if input_item.isdigit():
            analysis_id = "AN{}".format(input_item.zfill(6))
            yield GenericMWURL({
                "context": "study",
                "input_item": "analysis_id",
                "input_value": analysis_id,
                "output_item": "mwtab",
                "output_format": output_format
            }, base_url).url
        elif re.match(r'(AN[0-9]{6}$)', input_item):
            yield GenericMWURL({
                "context": "study",
                "input_item": "analysis_id",
                "input_value": input_item,
                "output_item": "mwtab",
                "output_format": output_format
            }, base_url).url
        elif re.match(r'(ST[0-9]{1,6}$)', input_item):
            yield GenericMWURL({
                "context": "study",
                "input_item": "study_id",
                "input_value": input_item,
                "output_item": "mwtab",
                "output_format": output_format
            }, base_url).url


def generate_urls(input_items, base_url=BASE_URL, **kwds):
    """
    Method for creating a generator which yields validated Metabolomics Workbench REST urls.

    :param list input_items: List of Metabolomics Workbench input values for mwTab files.
    :param str base_url: Base url to Metabolomics Workbench REST API.
    :param dict kwds: Keyword arguments of Metabolomics Workbench URL Path items.
    :return: Metabolomics Workbench REST URL string(s).
    :rtype: :py:class:`str`
    """
    for input_item in input_items:
        params = dict(kwds)
        params["input_item"] = input_item
        yield GenericMWURL(params, base_url).url


class GenericMWURL(object):
    """GenericMWURL class that stores and validates parameters specifying a Metabolomics Workbench REST URL.

    Metabolomics REST API requests are performed using URL requests in the form of
        https://www.metabolomicsworkbench.org/rest/context/input_specification/output_specification

        where:
            if context = "study" | "compound" | "refmet" | "gene" | "protein"
                input_specification = input_item/input_value
                output_specification = output_item/[output_format]
            elif context = "moverz"
                input_specification = input_item/input_value1/input_value2/input_value3
                    input_item = "LIPIDS" | "MB" | "REFMET"
                    input_value1 = m/z_value
                    input_value2 = ion_type_value
                    input_value3 = m/z_tolerance_value
                output_specification = output_format
                    output_format = "txt"
            elif context =  "exactmass"
                input_specification = input_item/input_value1/input_value2
                    input_item = "LIPIDS" | "MB" | "REFMET"
                    input_value1 = LIPID_abbreviation
                    input_value2 = ion_type_value
                output_specification = None

    """
    context = {
        "study": {
            "input_item": {
                'study_id', 'study_title', 'institute', 'last_name', 'analysis_id', 'metabolite_id'
            },
            "output_item": {
                'summary', 'factors', 'analysis', 'metabolites', 'mwtab', 'source', 'species', 'disease',
                'number_of_metabolites', 'data', 'datatable', 'untarg_studies', 'untarg_factors', 'untarg_data'
            }
        },

        "compound": {
            "input_item": {
                'regno', 'formula', 'inchi_key', 'lm_id', 'pubchem_cid', 'hmdb_id', 'kegg_id', 'chebi_id', 'metacyc_id',
                'abbrev'
            },
            "output_item": {
                'all', 'regno', 'formula', 'exactmass', 'inchi_key', 'name', 'sys_name', 'smiles', 'lm_id',
                'pubchem_cid', 'hmdb_id', 'kegg_id', 'chebi_id', 'metacyc_id', 'classification', 'molfile', 'png'
            }
        },

        "refmet": {
            "input_item": {
                'all', 'match', 'name', 'inchi_key', 'regno', 'pubchem_cid', 'formula', 'main_class', 'sub_class'
            },
            "output_item": {
                'all', 'name', 'inchi_key', 'regno', 'pubchem_cid', 'exactmass', 'formula', 'synonyms', 'sys_name',
                'main_class', 'sub_class'
            }
        },

        "gene": {
            "input_item": {
                'mgp_id', 'gene_id', 'gene_name', 'gene_symbol', 'taxid'
            },
            "output_item": {
                'all', 'lmp_id', 'mgp_id', 'gene_name', 'gene_symbol', 'gene_synonyms', 'alt_names', 'chromosome',
                'map_location', 'summary', 'taxid', 'species', 'species_long'
            }
        },

        "protein": {
            "input_item": {
                'mgp_id', 'gene_id', 'gene_name', 'gene_symbol', 'taxid', 'mrna_id', 'refseq_id', 'protein_gi',
                'uniprot_id', 'protein_entry', 'protein_name'
            },
            "output_item": {
                'all', 'mgp_id', 'gene_id', 'gene_name', 'gene_symbol', 'taxid', 'species', 'species_long', 'mrna_id',
                'refseq_id', 'protein_gi', 'uniprot_id', 'protein_entry', 'protein_name', 'seqlength', 'seq',
                'is_identical_to'
            }
        },

        "moverz": {
            "input_item": {
                'LIPIDS', 'MB', 'REFMET'
            },
            "ion_type_value": {
                'M+H', 'M+H-H2O', 'M+2H', 'M+3H', 'M+4H', 'M+K', 'M+2K', 'M+Na', 'M+2Na', 'M+Li', 'M+2Li', 'M+NH4',
                'M+H+CH3CN', 'M+Na+CH3CN', 'M.NaFormate+H', 'M.NH4Formate+H', 'M.CH3', 'M.TMSi', 'M.tBuDMSi', 'M-H',
                'M-H-H2O', 'M+Na-2H', 'M+K-2H', 'M-2H', 'M-3H', 'M4H', 'M.Cl', 'M.F', 'M.HF2', 'M.OAc', 'M.Formate',
                'M.NaFormate-H', 'M.NH4Formate-H', 'Neutral'
            }
        },

        "exactmass": {
            "LIPID_abbreviation": {
                'ArthroCer', 'asialoGM2Cer', 'CAR', 'CE', 'Cer', 'CerP', 'CoA', 'DG', 'DGDG', 'FA', 'GalCer', 'GB3Cer',
                'GlcCer', 'GM3Cer', 'GM4Cer', 'iGB3Cer', 'LacCer', 'Lc3Cer', 'Manb1-4GlcCer', 'MG', 'MGDG', 'MolluCer',
                'PA', 'PC', 'PE', 'PE-Cer', 'PG', 'PGP', 'PI', 'PI-Cer', 'PIP', 'PIP2', 'PIP3', 'PS', 'SM', 'SQDG', 'TG'
            },
            "ion_type_value": {
                'Neutral', 'M+H', 'M+H-H2O', 'M+2H', 'M+3H', 'M+4H', 'M+K', 'M+2K', 'M+2K-H', 'M+Na', 'M+2Na',
                'M+2Na-H',
                'M+Li', 'M+2Li', 'M+Ag', 'M+NH4', 'M-H', 'M-CH3', 'M2H', 'M-3H', 'M-4H', 'M.Cl', 'M.OAc', 'M.Formate'
            }
        }
    }

    def __init__(self, rest_params, base_url=BASE_URL):
        """URL initializer.

        :param dict rest_params: Dictionary of Metabolomics Workbench URL Path items.
        :param str base_url: Base url to Metabolomics Workbench REST API.
        """
        self.rest_params = rest_params
        self.base_url = base_url
        self._validate()
        self.url = self._create_url()

    def _validate(self):
        """Validate URL parameters. Raises error if self.rest_params is lacking a "context" keyword. Sub-functions raise
        error if missing or invalid parameter(s) in self.rest_params.

        :return: None
        :rtype: :py:obj:`None`
        """
        if not self.rest_params["context"] in self.context.keys():
            raise KeyError("Error: Invalid/missing context")
        elif self.rest_params["context"] in {"study", "compound", "refmet", "gene", "protein"}:
            self._validate_generic()
        elif self.rest_params["context"] == "moverz":
            self._validate_moverz()
        elif self.rest_params["context"] == "exactmass":
            self._validate_exactmass()

    def _validate_generic(self):
        """Validates keyword arguments for study, compound, refmet, gene, and protein contexts. Raises error if missing
        or invalid parameter(s) in self.rest_params.

        context = "study"
        input_item = "study_id" | "study_title" | "institute" | "last_name" | "analysis_id" | "metabolite_id"
        input_value = input_value
        output_item = "summary" | "factors" | "analysis" | "metabolites" | "mwtab" | "source" | "species" | "disease" |
                      "number_of_metabolites" | "data" | "datatable" | "untarg_studies" | "untarg_factors" |
                      "untarg_data"
        output_format = "txt" | "json" [Default: "json"]

        context = "compound"
        input_item = "regno" | "formula" | "inchi_key" | "lm_id" | "pubchem_cid" | "hmdb_id" | "kegg_id" | "chebi_id" |
                      "metacyc_id" | "abbrev"
        input_value = input_value
        output_item = "all" | "regno" | "formula" | "exactmass" | "inchi_key" | "name" | "sys_name" | "smiles" |
                      "lm_id" | "pubchem_cid" | "hmdb_id" | "kegg_id" | "chebi_id" | "metacyc_id" | "classification" |
                      "molfile" | "png" | "regno,formula,exactmass,..."
        output_format = "txt" | "json" [Default: "json"]

        context = "refmet"
        input_item = "all" | "match" | "name" | "inchi_key" | "regno" | "pubchem_cid" | "formula" | "main_class" |
                     "sub_class"
        input_value = input_value
        output_item = "all" | "name" | "inchi_key" | "regno" | "pubchem_cid" | "exactmass" | "formula" | "synonyms" |
                       "sys_name" | "main_class" | "sub_class" | "name,inchi_key,regno,..."
        output_format = "txt" | "json" [Default: "json"]

        context = "gene"
        input_item = "mgp_id" | "gene_id" | "gene_name" | "gene_symbol" | "taxid"
        input_value = input_value
        output_item = "all" | "lmp_id" | "mgp_id" | "gene_name" | "gene_symbol" | "gene_synonyms" | "alt_names" |
                      "chromosome" | "map_location" | "summary" | "taxid" | "species" | "species_long" |
                      "mgp_id,gene_id,gene_name,..."
        output_format = "txt" | "json" [Default: "json"]

        context = "protein"
        input_item = "mgp_id" | "gene_id" | "gene_name" | "gene_symbol" | "taxid" | "mrna_id" | "refseq_id" |
                     "protein_gi" | "uniprot_id" | "protein_entry" | "protein_name"
        input_value = input_value
        output_item = "all" | "mgp_id" | "gene_id" | "gene_name" | "gene_symbol" | "taxid" | "species" |
                      "species_long" | "mrna_id" | "refseq_id" | "protein_gi" | "uniprot_id" | "protein_entry" |
                      "protein_name" | "seqlength" | "seq" | "is_identical_to" | "mgp_id,gene_id,gene_name,..."
        output_format = "txt" | "json" [Default: "json"]

        Uses static method self._validate_id() to validate input_value parameter.

        :return: None
        :rtype: :py:obj:`None`
        """
        keywords = {"input_item", "input_value", "output_item"}
        # validate all required parameters are present
        if not all(k in self.rest_params.keys() for k in keywords):
            raise KeyError("Missing input item(s): " + str(keywords.difference(self.rest_params.keys())))

        # validate input_item
        if not any(k in self.rest_params["input_item"] for k in self.context[self.rest_params["context"]]["input_item"]):
            raise ValueError("Invalid input item")

        # validate input_item
        self._validate_input(self.rest_params["input_item"], self.rest_params["input_value"])

        # validate output_item(s)
        if type(self.rest_params["output_item"]) == list:
            if self.rest_params["context"] == "study":
                raise ValueError("Invalid output items. Study only takes a single output item.")
            elif not all(k in self.context[self.rest_params["context"]]["output_item"] for k in self.rest_params["output_item"]):
                raise ValueError("Invalid output item(s): " +
                                 str(set(self.rest_params["output_item"]).difference(self.context[self.rest_params["context"]]["output_item"])))
        elif not any(k in self.rest_params["output_item"] for k in self.context[self.rest_params["context"]]["output_item"]):
            raise ValueError("Invalid output item")

    def _validate_moverz(self):
        """Validate keyword arguments for moverz context. Raises error if missing or invalid parameter(s) in
        self.rest_params.

        context = "moverz"
        input_item = "LIPIDS" | "MB" | "REFMET"
        input_value1 = m/z_value
        input_value2 = ion_type_value
        input_value3 = m/z_tolerance_value
        output_format = "txt"

        m/z_value range: 50-2000
        See CONTEXT variable for supported ion type values.
        m/z_tolerance_value range: 0.0001-1

        :return: None
        :rtype: :py:obj:`None`
        """
        keywords = {"input_item", "m/z_value", "ion_type_value", "m/z_tolerance_value"}
        if not all(k in self.rest_params.keys() for k in keywords):
            raise KeyError("Missing input item(s): " + str(keywords.difference(self.rest_params.keys())))
        elif not any(k in self.rest_params["input_item"] for k in self.context["moverz"]["input_item"]):
            raise ValueError("Invalid input item")
        elif not 50 <= float(self.rest_params["m/z_value"]) <= 2000:
            raise ValueError("m/z value outside of range: 50-2000")
        elif not self.rest_params["ion_type_value"] in self.context["moverz"]["ion_type_value"]:
            raise ValueError("Invalid ion type value")
        elif not 0.0001 <= float(self.rest_params["m/z_tolerance_value"]) <= 1:
            raise ValueError("m/z tolerance value outside of range: 0.0001-1")

    def _validate_exactmass(self):
        """Validate keyword arguments for exactmass context. Raises error if missing or invalid parameter(s) in
        self.rest_params.

        context = "exactmass"
        LIPID_abbreviation = ...
        ion_type_value = ...

        See :class:`~mwtab.mwrest.GenericMWURL.context` variable for full list of possible values for LIPID_abbreviation
        and ion_type_value.

        :return: None
        :rtype: :py:obj:`None`
        """
        keywords = {"LIPID_abbreviation", "ion_type_value"}
        if not all(k in self.rest_params.keys() for k in keywords):
            raise KeyError("Missing input item(s): " + str(keywords.difference(self.rest_params.keys())))
        elif not any(k in self.rest_params["LIPID_abbreviation"] for k in self.context["exactmass"]["LIPID_abbreviation"]):
            raise ValueError("Invalid LIPID abbreviation")
        elif not self.rest_params['ion_type_value'] in self.context["exactmass"]["ion_type_value"]:
            raise ValueError("Invalid ion type value")

    @staticmethod
    def _validate_input(input_item, input_value):
        """Validate keyword arguments for input item where an id is used (ie. study, compound, refmet, gene, and
        protein). If invalid, raises value error.

        :param str input_item: String representing the input item from the input specification.
        :param str input_value: String representing the input vlaue from the input specification.
        :return: None
        :rtype: :py:obj:`None`
        """
        if input_item == "study_id":
            # allows for pulling a range of entries (ie. ST0001 pulls studies 100-199)
            if not re.match(r"(ST[0-9]{0,6}$)", input_value):
                raise ValueError("Invalid Metabolomics Workbench (MW) study ID (ST<6-digit integer>)")
        elif input_item in ["study_title", "institute", "last_name",
                            "formula",
                            "gene_name", "gene_symbol", "protein_entry", "protein_name"]:
            if not type(input_value) == str:
                raise ValueError("Invalid {} (<string>)".format(input_item.replace("_", " ")))
        elif input_item == "analysis_id":
            if not re.match(r'(AN[0-9]{6}$)', input_value):
                raise ValueError("Invalid Metabolomics Workbench analysis ID for a study (AN<6-digit integer>)")
        elif input_item == "metabolite_id":
            if not re.match(r"(ME[0-9]{6}$)", input_value):
                raise ValueError("Invalid Metabolomics Workbench metabolite ID for a study (ME<6-digit integer>)")

        elif input_item == "regno":
            if not input_value.isdigit():
                raise ValueError("Invalid Metabolomics Workbench Metabolite database ID (<integer>)")
        elif input_item == "inchi_key":
            if not re.match(r"([A-Z,\-]{27}$)", input_value):
                raise ValueError("Invalid InChIKey (27-character string)")
        elif input_item == "lm_id":
            if not re.match(r"(LM[A-Z]{2}[0-9]{8,10}$)", input_value):
                raise ValueError("Invalid LIPID MAPS ID (LM<2-character LIPID MAPS category><8-10 character string>)")
        elif input_item == "pubchem_cid":
            if not input_value.isdigit():
                raise ValueError("Invalid PubChem Compound ID (<integer>)")
        elif input_item == "hmdb_id":
            if not re.match(r"(HMDB[0-9]+$)", input_value):
                raise ValueError("Invalid Human Metabolome Database ID (HMDB<integer>)")
        elif input_item == "kegg_id":
            if not re.match(r"(CO[0-9]+$)", input_value):
                raise ValueError("Invalid KEGG compound ID (CO<integer>)")
        elif input_item == "chebi_id":
            if not input_value.isdigit():
                raise ValueError("Invalid ChEBI compound id (<integer>)")
        # TODO: update the following two input types.
        elif input_item == "metacyc_id":
            if not type(input_value) == str:
                raise ValueError("Invalid METACYC compound ID (<string>)")
        elif input_item == "abbrev":
            if not type(input_value) == str:
                raise ValueError("Invalid : Lipid bulk abbreviation (<string>)")

        elif input_item == "mgp_id":
            if not re.match(r"(MGP[0-9]{6}$)", input_value):
                raise ValueError("Invalid Human Metabolome Gene/Protein (MGP) database gene ID (MGP<6-digit integer>)")
        elif input_item == "gene_id":
            if not input_value.isdigit():
                raise ValueError("Invalid Entrez gene ID (<integer>)")
        elif input_item == "taxid":
            if not input_value.isdigit():
                raise ValueError("Invalid NCBI taxonomy ID (<integer>)")
        elif input_item == "mrna_id":
            if not re.match(r"(NM_[0-9]+$)", input_value):
                raise ValueError("Invalid mRNA ID (NM_<integer>)")
        elif input_item == "refseq_id":
            if not re.match(r"(NP_[0-9]+$)", input_value):
                raise ValueError("Invalid mRNA ID (NP_<integer>)")
        elif input_item == "protein_gi":
            if not input_value.isdigit():
                raise ValueError("Invalid NCBI protein GI (<integer>)")
        elif input_item == "uniprot_id":
            # regex from https://www.uniprot.org/help/accession_numbers
            if not re.match(r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}", input_value):
                raise ValueError("Invalid UniProt ID (see uniprot.org/help/accession_numbers)")

    def _create_url(self):
        """Method for constructing a formatted Metabolomics Workbench REST API URL from the given class parameters.

        :return: Formatted Metabolomics Workbench REST API URL.
        :rtype: :py:class:`str`
        """
        if self.rest_params["context"] in {"study", "compound", "refmet", "gene", "protein"}:
            # allows for url to include or not include output_format
            if self.rest_params.get("output_format"):
                return self.base_url + "/".join([self.rest_params.get(p) for p in [
                    "context", "input_item", "input_value", "output_item", "output_format"
                ]])
            else:
                return self.base_url + "/".join([self.rest_params.get(p) for p in [
                    "context", "input_item", "input_value", "output_item"
                ]])

        elif self.rest_params["context"] == "moverz":
            rest_params = [self.rest_params.get(p) for p in [
                "context", "input_item", "m/z_value", "ion_type_value", "m/z_tolerance_value"
            ]]
            rest_params.append("txt")
            return self.base_url + "/".join(rest_params)

        elif self.rest_params["context"] == "exactmass":
            return self.base_url + "/".join([self.rest_params.get(p) for p in [
                "context", "LIPID_abbreviation", "ion_type_value"
            ]])


class MWRESTFile(object):
    """MWRESTFile class that stores data from a single file download through Metabolomics Workbench's REST API.

    Mirrors :class:`~mwtab.mwtab.MWTabFile`.
    """

    def __init__(self, source):
        """File initializer.

        :param str source: Source a `MWRESTFile` instance was created from.
        """
        self.source = source
        self.text = ""

    def read(self, filehandle):
        """Read data into a :class:`~mwtab.mwrest.MWRESTFile` instance.

        :param filehandle: file-like object.
        :type filehandle: :py:class:`io.TextIOWrapper`, :py:class:`gzip.GzipFile`,
                          :py:class:`bz2.BZ2File`, :py:class:`zipfile.ZipFile`
        :return: None
        :rtype: :py:obj:`None`
        """
        self.text = filehandle.read().decode("utf-8")
        # input_str = input_str.replace("\r\n", "\n")
        # self.text = re.sub(r"</br>", "", self.text)  # included to remove remaining HTML tags
        filehandle.close()

    def write(self, filehandle):
        """Write :class:`~mwtab.mwrest.MWRESTFile` data into file.

        :param filehandle: file-like object.
        :type filehandle: :py:class:`io.TextIOWrapper`
        :return: None
        :rtype: :py:obj:`None`
        """
        try:
            filehandle.write(self.text)
        except IOError:
            raise IOError('"filehandle" parameter must be writable.')
        filehandle.close()

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
