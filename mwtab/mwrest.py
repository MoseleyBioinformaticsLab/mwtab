#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
mwtab.mwrest
~~~~~~~~~~~~

This module provides routines for accessing the ``mwTab`` REST API.
"""

from collections import OrderedDict

MWREST = "https://www.metabolomicsworkbench.org/rest/{}/"
CONTEXT = {
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
    },

    'exactmass': {}
}


class GenericMWURL(OrderedDict):

    def __init__(self, **kwds):
        super(GenericMWURL, self).__init__(**kwds)
        self._validate()

    def _validate(self):
        if not all(k in self.keys() for k in ['context', 'input item']):
            print("Error: Missing context or input item")
        elif not self['context'] in CONTEXT.keys():
            print("Error: Invalid input item")
        elif not self['input item'] in CONTEXT[self['context']]:
            print("Error: Invalid input item")
        pass
