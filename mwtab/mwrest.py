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


def analysis_ids():
    """
    Method for generating a list of urls for every current analysis in Metabolomics Workbench.

    :return: Urls to every Metabolomics Workbench analysis.
    :rtype: :py:class:`str`
    """
    st_an_dict = _pull_study_analysis()
    analyses = list()
    [analyses.extend(st_an_dict[k]) for k in st_an_dict.keys()]

    if VERBOSE:
        print("Found {} analysis file to download.".format(len(analyses)))

    return analyses


def study_ids():
    """
    Method for generating a list of urls for every current study in Metabolomics Workbench.

    :return: Urls to every Metabolomics Workbench study.
    :rtype: :py:class:`str`
    """
    st_an_dict = _pull_study_analysis()
    studies = list(st_an_dict.keys())

    if VERBOSE:
        print("Found {} study file to download.".format(len(studies)))

    return studies


def _pull_study_analysis():
    """
    Method for requesting a JSON string containing all study ids and analysis ids from Metabolomics Workbench's REST
    API. Requests a JSON file which contains a list of studies and their accompanying analyses. The JSON file is
    converted into a python object (dict) which can then be parsed to create a dictionary wit hthe form study id (key):
    analysis id(s) (values).

    :return: Dictionary of study ids (keys) and lists of analyses (value).
    :rtype: :py:class:`dict`
    """
    url = GenericMWURL(
        **{'context': 'study', 'input item': 'study_id', 'input value': 'ST', 'output item': 'analysis'}
    ).url
    mwrestfile = next(fileio.read_mwrest(url, **{'convertJSON': True}))
    json_object = mwrestfile._is_json(mwrestfile.text)

    study_analysis_dict = dict()
    for k in json_object.keys():
        if study_analysis_dict.get(json_object[k]['study_id']):
            study_analysis_dict[json_object[k]['study_id']].append(json_object[k]['analysis_id'])
        else:
            study_analysis_dict[json_object[k]['study_id']] = [json_object[k]['analysis_id']]

    return study_analysis_dict


def generate_mwtab_urls(input_items, output_format='txt'):
    """
    Method for generating URLS to be used to retrieve `mwtab` files for analyses and
    studies through the REST API of the Metabolomics Workbench database.

    :param list input_items:
    :param str output_format:
    :return: Metabolomics Workbench REST URL string.
    :rtype: :py:class:`str`
    """
    for input_item in input_items:
        if input_item.isdigit():
            analysis_id = "AN{}".format(input_item.zfill(6))
            yield GenericMWURL(**{
                'context': 'study',
                'input item': 'analysis_id',
                'input value': analysis_id,
                'output item': 'mwtab',
                'output format': output_format
            }).url
        elif re.match(r'(AN[0-9]{6}$)', input_item):
            yield GenericMWURL(**{
                'context': 'study',
                'input item': 'analysis_id',
                'input value': input_item,
                'output item': 'mwtab',
                'output format': output_format
            }).url
        elif re.match(r'(ST[0-9]{1,6}$)', input_item):
            yield GenericMWURL(**{
                'context': 'study',
                'input item': 'study_id',
                'input value': input_item,
                'output item': 'mwtab',
                'output format': output_format
            }).url


def generate_urls(input_items, **kwds):
    """
    Method for creating a generator which yields validated Metabolomics Workbench REST urls.

    :param list input_items:
    :param dict kwds:
    :return: Metabolomics Workbench REST URL string.
    :rtype: :py:class:`str`
    """
    for input_item in input_items:
        params = dict(kwds)
        params["input_item"] = input_item
        yield GenericMWURL(**params).url


class GenericMWURL(OrderedDict):
    """GenericMWURL class that stores Metabolomics Workbench REST API keyword argument data in the form of
    :py:class:`collections.OrderedDict`.

    Metabolomics REST API requests are performed using URL requests in the form of
        https://www.metabolomicsworkbench.org/rest/<context>/<input specification>/<output specification>

        where:
            <context> = study | compound | refmet | gene | protein | moverz | exactmass
            <input specification> = <input item>/<input value>
            <output specification> = <output item>/[<output format>]
    """
    base_mwrest_url = "https://www.metabolomicsworkbench.org/rest/"
    context = {
        'study': {
            'input item': {
                'study_id', 'study_title', 'institute', 'last_name', 'analysis_id', 'metabolite_id'
            },
            'output item': {
                'summary', 'factors', 'analysis', 'metabolites', 'mwtab', 'source', 'species', 'disease',
                'number_of_metabolites', 'data', 'datatable', 'untarg_studies', 'untarg_factors', 'untarg_data'
            }
        },

        'compound': {
            'input item': {
                'regno', 'formula', 'inchi_key', 'lm_id', 'pubchem_cid', 'hmdb_id', 'kegg_id', 'chebi_id', 'metacyc_id',
                'abbrev'
            },
            'output item': {
                'all', 'regno', 'formula', 'exactmass', 'inchi_key', 'name', 'sys_name', 'smiles', 'lm_id',
                'pubchem_cid', 'hmdb_id', 'kegg_id', 'chebi_id', 'metacyc_id', 'classification', 'molfile', 'png'
            }
        },

        'refmet': {
            'input item': {
                'all', 'match', 'name', 'inchi_key', 'regno', 'pubchem_cid', 'formula', 'main_class', 'sub_class'
            },
            'output item': {
                'all', 'name', 'inchi_key', 'regno', 'pubchem_cid', 'exactmass', 'formula', 'synonyms', 'sys_name',
                'main_class', 'sub_class'
            }
        },

        'gene': {
            'input item': {
                'mgp_id', 'gene_id', 'gene_name', 'gene_symbol', 'taxid'
            },
            'output item': {
                'all', 'lmp_id', 'mgp_id', 'gene_name', 'gene_symbol', 'gene_synonyms', 'alt_names', 'chromosome',
                'map_location', 'summary', 'taxid', 'species', 'species_long'
            }
        },

        'protein': {
            'input item': {
                'mgp_id', 'gene_id', 'gene_name', 'gene_symbol', 'taxid', 'mrna_id', 'refseq_id', 'protein_gi',
                'uniprot_id', 'protein_entry', 'protein_name'
            },
            'output item': {
                'all', 'mgp_id', 'gene_id', 'gene_name', 'gene_symbol', 'taxid', 'species', 'species_long', 'mrna_id',
                'refseq_id', 'protein_gi', 'uniprot_id', 'protein_entry', 'protein_name', 'seqlength', 'seq',
                'is_identical_to'
            }
        },

        'moverz': {
            'input item': {
                'LIPIDS', 'MB', 'REFMET'
            },
            'ion type value': {
                'M+H', 'M+H-H2O', 'M+2H', 'M+3H', 'M+4H', 'M+K', 'M+2K', 'M+Na', 'M+2Na', 'M+Li', 'M+2Li', 'M+NH4',
                'M+H+CH3CN', 'M+Na+CH3CN', 'M.NaFormate+H', 'M.NH4Formate+H', 'M.CH3', 'M.TMSi', 'M.tBuDMSi', 'M-H',
                'M-H-H2O', 'M+Na-2H', 'M+K-2H', 'M-2H', 'M-3H', 'M4H', 'M.Cl', 'M.F', 'M.HF2', 'M.OAc', 'M.Formate',
                'M.NaFormate-H', 'M.NH4Formate-H', 'Neutral'
            }
        },

        'exactmass': {
            'LIPID abbreviation': {
                'ArthroCer', 'asialoGM2Cer', 'CAR', 'CE', 'Cer', 'CerP', 'CoA', 'DG', 'DGDG', 'FA', 'GalCer', 'GB3Cer',
                'GlcCer', 'GM3Cer', 'GM4Cer', 'iGB3Cer', 'LacCer', 'Lc3Cer', 'Manb1-4GlcCer', 'MG', 'MGDG', 'MolluCer',
                'PA', 'PC', 'PE', 'PE-Cer', 'PG', 'PGP', 'PI', 'PI-Cer', 'PIP', 'PIP2', 'PIP3', 'PS', 'SM', 'SQDG', 'TG'
            },
            'ion type value': {
                'Neutral', 'M+H', 'M+H-H2O', 'M+2H', 'M+3H', 'M+4H', 'M+K', 'M+2K', 'M+2K-H', 'M+Na', 'M+2Na',
                'M+2Na-H',
                'M+Li', 'M+2Li', 'M+Ag', 'M+NH4', 'M-H', 'M-CH3', 'M2H', 'M-3H', 'M-4H', 'M.Cl', 'M.OAc', 'M.Formate'
            }
        }
    }

    def __init__(self, **kwds):
        """File initializer.

        :param dict kwargs: Dictionary of Metabolomics Workbench URL Path items.
        """
        super(GenericMWURL, self).__init__(**kwds)
        self.base_url = kwds.get("base url") or self.base_mwrest_url
        self.url = self._validate()

    def _validate(self):
        """Validate keyword arguments.

        :return: URL string.
        :rtype: :py:class:`str`
        """
        if not self['context'] in self.context.keys():
            raise KeyError("Error: Invalid/missing context")
        elif self['context'] in {'study', 'compound', 'refmet', 'gene', 'protein'}:
            return self._validate_generic()
        elif self['context'] == 'moverz':
            return self._validate_moverz()
        elif self['context'] == 'exactmass':
            return self._validate_exactmass()

    def _validate_generic(self):
        """Validate keyword arguments for study, compound, refmet, gene, and protein context.
        If valid, generates REST URL.

        <context> = study
        <input item> = study_id | study_title | institute | last_name | analysis_id | metabolite_id
        <input value> = <input item value>
        <output item> = summary | factors | analysis | metabolites | mwtab | source | species | disease |
            number_of_metabolites | data | datatable | untarg_studies | untarg_factors | untarg_data
        <output format> = txt | json [Default: json]

        <context> = compound
        <input item> = regno | formula | inchi_key | lm_id | pubchem_cid | hmdb_id | kegg_id | chebi_id |
            metacyc_id | abbrev
        <input value> = <input item value>
        <output item> = all | regno | formula | exactmass | inchi_key | name | sys_name | smiles | lm_id |
            pubchem_cid | hmdb_id | kegg_id | chebi_id | metacyc_id | classification | molfile | png |
            regno,formula,exactmass,...
        <output format> = txt | json [Default: json]

        <context> = refmet
        <input item> = all | match | name | inchi_key | regno | pubchem_cid | formula | main_class | sub_class
        <input value> = <input item value>
        <output item> = all | name | inchi_key | regno | pubchem_cid | exactmass | formula | synonyms | sys_name |
            main_class | sub_class | name,inchi_key,regno,...
        <output format> = txt | json

        <context> = gene
        <input item> = mgp_id | gene_id | gene_name | gene_symbol | taxid
        <input value> = <input item value>
        <output item> = all | lmp_id | mgp_id | gene_name | gene_symbol | gene_synonyms | alt_names | chromosome |
            map_location | summary | taxid | species | species_long | mgp_id,gene_id,gene_name,...
        <output format> = txt | json

        <context> = protein
        <input item> = mgp_id | gene_id | gene_name | gene_symbol | taxid | mrna_id | refseq_id | protein_gi |
            uniprot_id | protein_entry | protein_name
        <input value> = <input item value>
        <output item> = all | mgp_id | gene_id | gene_name | gene_symbol | taxid | species | species_long | mrna_id |
            refseq_id | protein_gi | uniprot_id | protein_entry | protein_name | seqlength | seq | is_identical_to |
            mgp_id,gene_id,gene_name,...
        <output format> = txt | json

        Uses static method self._validate_id() to validate input value.

        :return: URL string.
        :rtype: :py:class:`str`
        """
        keywords = {'input item', 'input value', 'output item'}
        if not all(k in self.keys() for k in keywords):
            raise KeyError("Missing input item(s): " + str(keywords.difference(self.keys())))
        elif not any(k in self['input item'] for k in self.context[self['context']]['input item']):
            raise ValueError("Invalid input item")
        elif self._validate_input(self['input item'], self['input value']):
            exit()
        elif type(self['output item']) == list:
            if self['context'] == 'study':
                raise ValueError("Invalid output items. Study only takes a single output item.")
            elif not all(k in self.context[self['context']]['output item'] for k in self['output item']):
                raise ValueError("Invalid output item(s): " +
                                 str(set(self['output item']).difference(self.context[self['context']]['output item'])))
            else:
                return self.base_url + '/'.join([
                    self.get('context'),
                    self.get('input item'),
                    str(self.get('input value')),
                    ','.join(self.get('output item')),
                    self.get('output format') or ''
                ])
        elif not any(k in self['output item'] for k in self.context[self['context']]['output item']):
            raise ValueError("Invalid output item")
        else:
            return self.base_url + '/'.join([self.get('context'), self.get('input item'), str(self.get('input value')),
                                      self.get('output item'), self.get('output format') or ''])

    def _validate_moverz(self):
        """Validate keyword arguments for moverz context. If valid, generates REST URL.

        <context> = moverz
        <input item> = LIPIDS | MB | REFMET
        <input value1> = <m/z value>
        <input value2> = <ion type value>
        <input value3> = <m/z tolerance value>
        <output format> = txt

        <m/z value> range: 50-2000
        See CONTEXT variable for supported ion type values.
        <m/z tolerance value> range: 0.0001-1

        :return: URL string.
        :rtype: :py:class:`str`
        """
        keywords = {'input item', 'm/z value', 'ion type value', 'm/z tolerance value'}
        if not all(k in self.keys() for k in keywords):
            raise KeyError("Missing input item(s): " + str(keywords.difference(self.keys())))
        elif not any(k in self['input item'] for k in self.context['moverz']['input item']):
            raise ValueError("Invalid input item")
        elif not 50 <= float(self['m/z value']) <= 2000:
            raise ValueError("m/z value outside of range: 50-2000")
        elif not self['ion type value'] in self.context['moverz']['ion type value']:
            raise ValueError("Invalid ion type value")
        elif not 0.0001 <= float(self['m/z tolerance value']) <= 1:
            raise ValueError("m/z tolerance value outside of range: 0.0001-1")
        else:
            # only supports txt output format
            return self.base_url + '/'.join([self['context'], self['input item'], str(self['m/z value']),
                                      self['ion type value'], str(self['m/z tolerance value']), 'txt'])

    def _validate_exactmass(self):
        """Validate keyword arguments for exactmass context. If valid, generates REST URL.

        <context> = exactmass
        <input value1> = <LIPID abbreviation>
        <input value2> = <ion type value>

        See CONTEXT variable for supported LIPID abbreviations and ion type values.

        :return: URL string.
        :rtype: :py:class:`str`
        """
        keywords = {'LIPID abbreviation', 'ion type value'}
        if not all(k in self.keys() for k in keywords):
            raise KeyError("Missing input item(s): " + str(keywords.difference(self.keys())))
        elif not any(k in self['LIPID abbreviation'] for k in self.context['exactmass']['LIPID abbreviation']):
            raise ValueError("Invalid LIPID abbreviation")
        elif not self['ion type value'] in self.context['exactmass']['ion type value']:
            raise ValueError("Invalid ion type value")
        else:
            # no required output format
            return self.base_url + '/'.join([self['context'], self['LIPID abbreviation'], self['ion type value']])

    @staticmethod
    def _validate_input(input_item, input_value):
        """Validate keyword arguments for input item where an id is used (ie. study, compound, refmet, gene, and
        protein). If invalid, raises value error.
        """
        if input_item == 'study_id':
            # allows for pulling a range of entries (ie. ST0001 pulls studies 100-199)
            if not re.match(r'(ST[0-9]{0,6}$)', input_value):
                raise ValueError("Invalid Metabolomics Workbench (MW) study ID (ST<6-digit integer>)")
        elif input_item in ['study_title', 'institute', 'last_name',
                            'formula',
                            'gene_name', 'gene_symbol', 'protein_entry', 'protein_name']:
            if not type(input_value) == str:
                raise ValueError("Invalid {} (<string>)".format(input_item.replace('_', ' ')))
        elif input_item == 'analysis_id':
            if not re.match(r'(AN[0-9]{6}$)', input_value):
                raise ValueError("Invalid Metabolomics Workbench analysis ID for a study (AN<6-digit integer>)")
        elif input_item == 'metabolite_id':
            if not re.match(r'(ME[0-9]{6}$)', input_value):
                raise ValueError("Invalid Metabolomics Workbench metabolite ID for a study (ME<6-digit integer>)")

        elif input_item == 'regno':
            if not input_value.isdigit():
                raise ValueError("Invalid Metabolomics Workbench Metabolite database ID (<integer>)")
        elif input_item == 'inchi_key':
            if not re.match(r'([A-Z,\-]{27}$)', input_value):
                raise ValueError("Invalid InChIKey (27-character string)")
        elif input_item == 'lm_id':
            if not re.match(r'(LM[A-Z]{2}[0-9]{8,10}$)', input_value):
                raise ValueError("Invalid LIPID MAPS ID (LM<2-character LIPID MAPS category><8-10 character string>)")
        elif input_item == 'pubchem_cid':
            if not input_value.isdigit():
                raise ValueError("Invalid PubChem Compound ID (<integer>)")
        elif input_item == 'hmdb_id':
            if not re.match(r'(HMDB[0-9]+$)', input_value):
                raise ValueError("Invalid Human Metabolome Database ID (HMDB<integer>)")
        elif input_item == 'kegg_id':
            if not re.match(r'(CO[0-9]+$)', input_value):
                raise ValueError("Invalid KEGG compound ID (CO<integer>)")
        elif input_item == 'chebi_id':
            if not input_value.isdigit():
                raise ValueError("Invalid ChEBI compound id (<integer>)")
        # TODO: update the following two input types.
        elif input_item == 'metacyc_id':
            if not type(input_value) == str:
                raise ValueError("Invalid METACYC compound ID (<string>)")
        elif input_item == 'abbrev':
            if not type(input_value) == str:
                raise ValueError("Invalid : Lipid bulk abbreviation (<string>)")

        elif input_item == 'mgp_id':
            if not re.match(r'(MGP[0-9]{6}$)', input_value):
                raise ValueError("Invalid Human Metabolome Gene/Protein (MGP) database gene ID (MGP<6-digit integer>)")
        elif input_item == 'gene_id':
            if not input_value.isdigit():
                raise ValueError("Invalid Entrez gene ID (<integer>)")
        elif input_item == 'taxid':
            if not input_value.isdigit():
                raise ValueError("Invalid NCBI taxonomy ID (<integer>)")
        elif input_item == 'mrna_id':
            if not re.match(r'(NM_[0-9]+$)', input_value):
                raise ValueError("Invalid mRNA ID (NM_<integer>)")
        elif input_item == 'refseq_id':
            if not re.match(r'(NP_[0-9]+$)', input_value):
                raise ValueError("Invalid mRNA ID (NP_<integer>)")
        elif input_item == 'protein_gi':
            if not input_value.isdigit():
                raise ValueError("Invalid NCBI protein GI (<integer>)")
        elif input_item == 'uniprot_id':
            # regex from https://www.uniprot.org/help/accession_numbers
            if not re.match(r'[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}', input_value):
                raise ValueError("Invalid UniProt ID (see uniprot.org/help/accession_numbers)")


class MWRESTFile(object):
    """MWRESTFile class that stores data from a single file download
    through Metabolomics Workbench's REST API.
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
        input_str = filehandle.read().decode('utf-8')
        self.text = input_str
        filehandle.close()

    def write(self, filehandle):
        """Write :class:`~mwtab.mwrest.MWRESTFile` data into file.

        :param filehandle: file-like object.
        :type filehandle: :py:class:`io.TextIOWrapper`
        :param str file_format: Format to use to write data.
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
