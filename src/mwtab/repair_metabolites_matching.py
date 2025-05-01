# -*- coding: utf-8 -*-
"""
Contains regular expressions specific to the repair command use case and the column_matching_attributes data structure.
"""

import pandas

from . import metabolites_regexes

make_list_regex = metabolites_regexes.make_list_regex



INTEGER = r'-?\d+'
FLOAT = r'-?\d*\.\d+'
SCIENTIFIC_NOTATION = r'-?\d*\.\d+E(-|\+)?\d+'
NUMS = r'(' + FLOAT + '|' + SCIENTIFIC_NOTATION + '|' + INTEGER + r')'
# LIST_OF_NUMS = r'(' + NUMS + r'\s*,\s*)+' + r'(' + NUMS + r'\s*|\s*)'
LIST_OF_NUMS = make_list_regex(NUMS, ',')
BRACKETED_LIST_OF_NUMS = r'\[' + LIST_OF_NUMS + r'\]'
PARENTHESIZED_LIST_OF_NUMS = r'\(' + LIST_OF_NUMS + r'\)'
LIST_OF_NUMS_UNDERSCORE = make_list_regex(NUMS, '_')
LIST_OF_NUMS_SLASH = make_list_regex(NUMS, '/')
POSITIVE_NUMS = NUMS.replace('-?', '')
POSITIVE_INTS = r'\d+'
LIST_OF_POS_INTS = make_list_regex(POSITIVE_INTS, ',')
LIST_OF_POS_INTS_OR = make_list_regex(POSITIVE_INTS, 'or')
LIST_OF_POS_INTS_BAR = make_list_regex(POSITIVE_INTS, '\|')
LIST_OF_POS_INTS_SLASH = make_list_regex(POSITIVE_INTS, '/')
LIST_OF_POS_INTS_SEMICOLON = make_list_regex(POSITIVE_INTS, ';')
LIST_OF_POS_INTS_SPACE = make_list_regex(POSITIVE_INTS, ' ')

POSITIVE_FLOATS = r'\d*.\d+'
POSITIVE_SCIENTIFIC_NOTATION = r'\d*\.\d*E(-|\+)?\d+'
POSITIVE_FLOAT_RANGE = POSITIVE_FLOATS + r'\s*(_|-)\s*' + POSITIVE_FLOATS
LIST_OF_POS_FLOATS_UNDERSCORE = make_list_regex(POSITIVE_FLOATS, '_')
POS_FLOAT_PAIRS = r'(' +  POSITIVE_FLOATS + r'(_|@)' + POSITIVE_FLOATS + r')'
LIST_OF_POS_FLOAT_PAIRS_UNDERSCORE = make_list_regex(POS_FLOAT_PAIRS, '_')
LIST_OF_POS_FLOAT_PAIRS_NO_SPACE = make_list_regex(POS_FLOAT_PAIRS, '')
LIST_OF_POS_FLOAT_PAIRS_MIXED = make_list_regex(POS_FLOAT_PAIRS, '(//|,)')
POS_INT_FLOAT_PAIR = r'(' +  POSITIVE_INTS + r'_' + POSITIVE_FLOATS + r')'

# ELEMENT_SYMBOL = r'[A-Z][a-z]?'
ELEMENT_SYMBOL = r'([BCFHIKNOPSUVWY]|[ISZ][nr]|[ACELP][ru]|A[cglmst]|B[aehikr]|C[adeflos]|D[bsy]|Es|F[elmr]|G[ade]|H[efgos]|Kr|L[aiv]|M[cdgnot]|N[abdehiop]|O[gs]|P[abdmot]|R[abe-hnu]|S[bcegim]|T[abcehilms]|Xe|Yb)'
ELEMENT_COUNT = r'([1-9]\d*)*'
FORMULA_ELEMENT = ELEMENT_SYMBOL + ELEMENT_COUNT
FORMULA = r'(' + FORMULA_ELEMENT + ')+'
LIST_OF_FORMULAS = make_list_regex(FORMULA, ',', True)
BRACKETED_LIST_OF_FORMULAS = r'\[' + LIST_OF_FORMULAS + r'\]'
# C12H22O11, C12 H22 O11, [13C]4H7NO4, [13]C6 H14 [15]N4 O2, [C13]C4H6O5, C21C(13)2H38N7O17P3S, 12C12+14N4+16O19+1H32
ISOTOPIC_NUM = r'\d+'
# Add dueterium.
ISOTOPIC_SYMBOL = ELEMENT_SYMBOL[0:-1] + r'|D)'
ISOTOPIC_ELEMENT = ISOTOPIC_SYMBOL + ELEMENT_COUNT
ISOTOPIC_FORMULA = r'(' + ISOTOPIC_ELEMENT + '|' + \
                   '\[' + ISOTOPIC_NUM + ISOTOPIC_SYMBOL + '\]' + ELEMENT_COUNT + '|' + \
                   '\[' + ISOTOPIC_NUM + '\]' + ISOTOPIC_ELEMENT + '|' + \
                   '\[' + ISOTOPIC_SYMBOL + ISOTOPIC_NUM + '\]' + ELEMENT_COUNT + '|' + \
                   ISOTOPIC_SYMBOL + '\(' + ISOTOPIC_NUM + '\)' + ELEMENT_COUNT + '|' + \
                   make_list_regex(ISOTOPIC_NUM + ISOTOPIC_SYMBOL + ELEMENT_COUNT, '\+') + ')+'
LIST_OF_ISOTOPIC_FORMULAS = make_list_regex(ISOTOPIC_FORMULA, ',', True)
BRACKETED_LIST_OF_ISOTOPIC_FORMULAS = r'\[' + LIST_OF_ISOTOPIC_FORMULAS + r'\]'
# C12H22O11+, [C12H22O11]+
CHARGE_FORMULA = r'(\[' + FORMULA + r'\](-|\+)' + r'|' + FORMULA + '(-|\+)' + r')'
# CH3(CH2)16COOH
GROUP_FORMULA = r'(' + FORMULA_ELEMENT + '|' + r'\(' + FORMULA + r'\)\d+' + r')+'

ORGANIC_ELEMENT_SYMBOL = r'([CHNOPS])'
ORGANIC_FORMULA_ELEMENT = ORGANIC_ELEMENT_SYMBOL + ELEMENT_COUNT
ORGANIC_FORMULA = r'(' + ORGANIC_FORMULA_ELEMENT + '){4,}'

SMILES_ELEMENT_SYMBOL = r'\d?[a-zA-Z][a-z]?'
# Note the , , , and  symbols at the end of the character set are '\x01', 'x02', '\x03', and '\x04'. 
# I don't think they are apart of the SMILES characters, but they appeared in the SMILES of some datasets.
# There are some SMILES in AN003143 like 'C[n]1c[n]cc1C[C@@H](NC(=O)CCN)C(O)=O |&1:7|'. 
# I don't think the end bit between '|' is part of a legit SMILES, but I added it to pass this dataset.
SMILES = r'(\[?' + SMILES_ELEMENT_SYMBOL + r'[0-9@+\-[\]()\\/%=#$.*]*)+' + r'( \|[0-9:&,w]+\|)?'
LIST_OF_SMILES_SEMICOLON = make_list_regex(SMILES, ';')

INCHIKEY = r'(InChIKey=)?[a-zA-Z]{14}-[a-zA-Z]{10}-[a-zA-Z]?'
INCHIKEY_OR_NULL = r'(' + INCHIKEY + r'|null|No record)'
LIST_OF_INCHIKEYS = r'(' + INCHIKEY + '\s*,\s*)+' + r'(' + INCHIKEY + r'\s*|\s*)'
LIST_OF_INCHIKEYS_SLASH = LIST_OF_INCHIKEYS.replace(',', '/')

# InChI=1 or InChI=1 where "1" is a version number and "S" means it's a standard InChI.
INCHI_HEAD = r'InChI='
INCHI_VERSION = r'\d+S?'
INCHI_FORMULA = r'/' + FORMULA
INCHI_SKELETAL_LAYER = r'/c(\d+([-,()])?)+'
# INCHI_HYDROGEN_LAYER = r'/h(\d+(-\d+)?,)*' + r'(\d+(-\d+)?H\d*,)*' + r'(\(H\d*(,\d+)+\))*'
# INCHI_HYDROGEN_LAYER = r'/h' + r'(' + make_list_regex(r'\d+(-\d+)?', ',') + r')?' + r'(' + make_list_regex(r'\d+(-\d+)?H\d*', ',') + r')?' + r'(\(H\d*(,\d+)+\))*'
INCHI_HYDROGEN_LAYER = r'/h' + r'(' + make_list_regex(r'\d+(-\d+)?H?\d*', ',') + r')?' + r'(\(H\d*-?(,\d+)+\))*'
INCHI_CHARGE_LAYER = r'/q(-|\+)\d+'
INCHI_PROTONATION_LAYER = r'/p(-|\+)\d+'
# INCHI_STEREOCHEMISTRY_LAYER = r'/t\d+(-|\+)(,\d+(-|\+))+'
INCHI_STEREOCHEMISTRY_LAYER = r'/t' + make_list_regex(r'(\d+(-|\+|\?|u)|M)', ',')
INCHI_STEREOCHEMISTRY_SUBLAYER1 = r'/m\d+'
INCHI_STEREOCHEMISTRY_SUBLAYER2 = r'/s\d+'
INCHI_FIXED_HYDROGEN_LAYER = r'/f(' + FORMULA + r')?' + \
                             r'(' + INCHI_HYDROGEN_LAYER + r')?' + \
                             r'(' + INCHI_STEREOCHEMISTRY_SUBLAYER2 + r')?' + \
                             r'(' + INCHI_CHARGE_LAYER + r')?'
INCHI_DOUBLEBOND_LAYER = r'/b' + make_list_regex(r'\d+(-|\+)\d+(-|\+)', ',')
# INCHI_ISOTOPIC_LAYER = r'/i\d+(-|\+)\d+(,\d+(-|\+)\d+)*' + r'(' + INCHI_STEREOCHEMISTRY_SUBLAYER2 + r')?'
# 1+1, 2H3, 2H, 1+1H3
INCHI_ISOTOPIC_LAYER = r'/i' + make_list_regex(r'(\d+(-|\+)\d+|\d+[A-Z]\d*|\d+(-|\+)\d+[A-Z]?\d*)', ',') + r'(' + INCHI_STEREOCHEMISTRY_SUBLAYER2 + r')?'

FULL_INCHI = r'(' + INCHI_HEAD + r')?' + \
             INCHI_VERSION + \
             r'(' + INCHI_FORMULA + r')?' + \
             r'(' + INCHI_SKELETAL_LAYER + r')?' + \
             r'(' + INCHI_HYDROGEN_LAYER + r')?' + \
             r'(' + INCHI_CHARGE_LAYER + r')?' + \
             r'(' + INCHI_PROTONATION_LAYER + r')?' + \
             r'(' + INCHI_DOUBLEBOND_LAYER + r')?' + \
             r'(' + INCHI_STEREOCHEMISTRY_LAYER + r')?' + \
             r'(' + INCHI_STEREOCHEMISTRY_SUBLAYER1 + r')?' + \
             r'(' + INCHI_STEREOCHEMISTRY_SUBLAYER2 + r')?' + \
             r'(' + INCHI_ISOTOPIC_LAYER + r')?' + \
             r'(' + INCHI_FIXED_HYDROGEN_LAYER + r')?'

CHEAP_INCHI = r'\s*' + r'(' + INCHI_HEAD + r')?' + INCHI_VERSION + r'((/c|/h|/i|/t|/m|/s|/q|/f|/p|/b|/[A-Z])(\S)+)+' + r'\s*'
LIST_OF_INCHI = make_list_regex(CHEAP_INCHI, ',', True)
BRACKETED_LIST_OF_INCHI = r'\[' + LIST_OF_INCHI + r'\]'

# M+H, [M+H]+, -H(-), M+AGN+H, M+Acid, (M-H)/2, [M+NH4] +_[M+Na]+, [M+H–H2O]+, [M-2H](2-)
# Cat, Hac, and Chol-head are wierd special cases.
ION_ELEMENTS = r'(\d?(m|M)' + r'|' + \
               r'(-|\+)?\d*' + FORMULA + r'\d*(\((-|\+)\))?' +  r'|' + \
               r'[ [a-zA-Z]*(A|a)cid' + r'|' + \
               r'Cat' + r'|' + r'Chol-head' + r'|' + r'Hac' + r'|' + r'\di' + r'|' + r'FA|NA|A' + r')'
ION = r'(' + make_list_regex(ION_ELEMENTS, '(-|\+)?') + r')'
ION = r'(' + ION + r'|' + \
      r'\[' + ION + r'\]\(?\d?\s?(-|\+)?\d?\)?' + r'|' + \
      ION + r'\](-|\+)?' + r'|' + \
      r'\(' + ION + r'\)((/\d)|(-|\+))?' + r')'
LIST_OF_IONS = make_list_regex(ION, ',')
LIST_OF_IONS_SPACE = make_list_regex(ION, ' ')
LIST_OF_IONS_UNDERSCORE = make_list_regex(ION, '_')
LIST_OF_IONS_MIXED = make_list_regex(ION, '(_| )')
LIST_OF_IONS_NO_DELIMITER = make_list_regex(ION, '')
BRACKETED_LIST_OF_IONS = r'\[' + make_list_regex(ION, ',', True) + r'\]'

# There are values like CA1511 that I can't confirm are KEGG values. You can't find them in the compound database, but they appear often.
# The prefixes CE, UP, Z, and U are the same. There is a common mistake of only having 4 numbers after C, so that is accounted for as well.
KEGG = r'(' + r'(cpd:)?[CDMGKRZU]0?\d{5}\?{0,2}' + r'|' + r'(DG|ko)\d{5}' + r'|' + r'(CA|CE|UP|C)\d{4}' + r'|' + r'(NA|n/a)' + r')'
LIST_OF_KEGG = make_list_regex(KEGG, ',')
LIST_OF_KEGG_SEMICOLON = make_list_regex(KEGG, ';')
LIST_OF_KEGG_SLASH = make_list_regex(KEGG, '/')
LIST_OF_KEGG_SPACE = make_list_regex(KEGG, ' ')
LIST_OF_KEGG_DOUBLE_SLASH = make_list_regex(KEGG, '//')
LIST_OF_KEGG_UNDERSCORE = make_list_regex(KEGG, '_')
LIST_OF_KEGG_HYPHEN = make_list_regex(KEGG, '-')
LIST_OF_KEGG_MIXED = make_list_regex(KEGG, '(/|,)')
LIST_OF_KEGG_BAR = make_list_regex(KEGG, '(\|)')

HMDB = r'(' + r'(HMDB|HDMB|YMDB|HMBD)\d+(\*|\?)?' + r'|' + r'n/a' + r')'
LIST_OF_HMDB = make_list_regex(HMDB, ',')
LIST_OF_HMDB_SLASH = make_list_regex(HMDB, '/')
LIST_OF_HMDB_AMPERSAND = make_list_regex(HMDB, '&')
LIST_OF_HMDB_SEMICOLON = make_list_regex(HMDB, ';')
LIST_OF_HMDB_SPACE = make_list_regex(HMDB, ' ')
LIST_OF_HMDB_UNDERSCORE = make_list_regex(HMDB, '_')
HMDB_INT = r'\d{,5}'
LIST_OF_HMDB_INTS = make_list_regex(HMDB_INT, ',')
LIST_OF_HMDB_INTS_SLASH = make_list_regex(HMDB_INT, '/')

LIPID_MAPS = r'(' + 'LM(PK|ST|GL|FA|SP|GP|PR|SL)[0-9A-Z]{8,10}\*?' + r'|' + r'(ST|FA|PR|GP|PK|GL|SP)\d{4,6}-' + FORMULA + r')'
LIST_OF_LMP = make_list_regex(LIPID_MAPS, ',')
LIST_OF_LMP_UNDERSCORE = make_list_regex(LIPID_MAPS, '_')
LIST_OF_LMP_SLASH = make_list_regex(LIPID_MAPS, '/')

# When using Excel, a CAS number can get mistaken for a date and it will automatically change the value.
DATE = r'\d{1,2}/\d{1,2}/(\d{4}|\d{2})'
# Note that only \d+-\d\d-\d is an actual CAS number, other things are there to catch common mistakes.
# CAS = r'(CAS: ?)?\d+-\d\d-0?\d' + r'|' + make_list_regex(r'\d+', '-') + r'|' + DATE
CAS = r'(CAS: ?)?\d+-\d\d-0?\d' + r'|' + DATE
LIST_OF_CAS = make_list_regex(CAS, ',')
LIST_OF_CAS_SEMICOLON = make_list_regex(CAS, ';')


class ColumnMatcher:
    def __init__(self, regex_search_strings: None|list[str] = None, 
                 not_regex_search_strings: None|list[str] = None, regex_search_sets: None|list[list[str]] = None,
                 in_strings: None|list[str] = None, not_in_strings: None|list[str] = None, 
                 in_string_sets: None|list[list[str]] = None, exact_strings: None|list[str] = None,
                 values_type: None|str = None, values_regex: None|str = None, inverse_values_regex: None|str = None):
        self.regex_search_strings = regex_search_strings if regex_search_strings else []
        self.not_regex_search_strings = not_regex_search_strings if not_regex_search_strings else []
        self.regex_search_sets = regex_search_sets if regex_search_sets else []
        self.in_strings = in_strings if in_strings else []
        self.not_in_strings = not_in_strings if not_in_strings else []
        self.in_string_sets = in_string_sets if in_string_sets else []
        self.exact_strings = exact_strings if exact_strings else []
        
        self.values_type = values_type if isinstance(values_type, str) else ''
        self.values_regex = values_regex if isinstance(values_regex, str) else ''
        self.inverse_values_regex = inverse_values_regex if isinstance(inverse_values_regex, str) else ''
    
    def series_values_match(self, series: pandas.Series) -> pandas.Series:
        """Return a mask for the series based on type and regex matching.
        
        "values_regex" and "inverse_values_regex" are mutually exclusive and "values_regex" will take precedence if both are given. 
        "values_type" and one of the regex parameters can both be used, the intermediate selectors are ANDed together.
        
        Args:
            series: series to match values based on type and/or regex.
        
        Returns:
            A pandas Series the same length as "series" with Boolean values that can be used to select the matching values in the series.
        """
        if self.values_regex:
            regex_match = series.str.fullmatch(self.values_regex, na=False)
        elif self.inverse_values_regex:
            regex_match = ~series.str.fullmatch(self.inverse_values_regex, na=True)
        else:
            regex_match = pandas.Series([True]*len(series), index=series.index)
        
        old_NAs = series.isna()
        column_to_numeric = pandas.to_numeric(series, errors='coerce')
        column_to_numeric_NAs = column_to_numeric.isna()
        new_NAs = column_to_numeric_NAs ^ old_NAs
        
        if self.values_type == "integer":
            # The top line will return True for values like '1.0', but the bottom line won't.
            # type_match = (column_to_numeric % 1 == 0) | ~new_NAs
            type_match = ~series.str.contains('.', regex=False, na=False) | old_NAs
        elif self.values_type == "numeric":
            type_match = ~new_NAs
        elif self.values_type == "non-numeric":
            type_match = new_NAs | old_NAs
        else:
            type_match = pandas.Series([True]*len(series), index=series.index)
        
        return regex_match & type_match
        


# TODO check with the Metabolomics Workbench and see if they have preferred harmonized names, not sure why Christian did 'moverz' and 'moverz_quant' separately.
# There might need to be 2 sets of dicts. One that is tighter detection on the names for name harmonization, 
# and another that is looser for validating values. For instance, not all column names currently detected for m/z should 
# be harmonized to moverz or whatever.
# Probably going to have to be aware of 2 column names matching at the same time, for example CLASS and CLASS_2 in AN004528.

# AN001956 has some values like ", HMDB00510" where the first value in the list isn't known. Need to decide if that's the proper way 
# that should be. Currently the regexes don't recognize this as a list, but maybe it should.

# Having a standard set of columns and trying to fill them in is a slightly different approach than 
# trying to classify and repair/harmonize known columns that already exist. Unfortuneately, I think 
# we have to do both. They require slightly different data structures, but a lof of the info is shared.
# I think I need to make the one for the classification and then the standard set one will be a subset/collection of those.
# For instance, the "identifier" column would not be a standard column, but could be used to fill in the "moverz" standard column.
column_matching_attributes = {    
    # 'MW structure' and 'MW regno' are IDs to the Metabolomics Workbench structure database.
    'moverz_quant': {
        # Other checks, non-negative, not all zeros
        # Need to look for the columns with values like mz_rt.
        'names': ['m/z', 'characteristic m/z', 'moverz quant', 'mz', 'quantified m/z'],
        'regex_search_strings': ['m/z', 'mz', 'moverz',  'mx'],
        'regex_search_sets': [],
        'not_regex_search_strings': ['id'],
        'in_strings': ['m.z', 'calcmz', 'medmz', 'm_z', 'obsmz', 'mass to charge', 'mass over z'],
        'in_string_sets': [],
        'not_in_strings': ['spec', 'pectrum', 'structure', 'regno', 'retention'],
        'exact_strings': [],
        'type': None,
        'regex': NUMS + r'|' + \
                 LIST_OF_NUMS + r'|' + \
                 LIST_OF_NUMS_UNDERSCORE + r'|' + \
                 NUMS + r'\s*/\s*' + NUMS + r'|' + \
                 r'(' + NUMS + r'\s*\(' + NUMS + r'\)' + '\s*;\s*)+' + r'(' + NUMS + r'\s*\(' + NUMS + r'\)' + r'\s*|\s*)' + r'|' + \
                 NUMS + r'(\s*>\s*|\s*<\s*)' + NUMS + r'|' + \
                 r'(' + NUMS + r'\s*)+' + r'|' + \
                 NUMS + r'\s*\(\s*' + NUMS + r'\s*\)' + r'|' + \
                 NUMS + r'\s*-\s*' + NUMS
        },
    'mass': {
        # Other checks, non-negative, not all zeros
        # Need to look for the columns with values like mz_rt.
        'names': ['mass', 'exact mass', 'monoisotopic mass'],
        'regex_search_strings': ['mass', 'quantmass', 'masses', 'mw', 'weight'],
        'regex_search_sets': [],
        'not_regex_search_strings': ['id'],
        'in_strings': ['exactmass', 'obsmass', 'calcmass', 'monoisotopicmass', 'molwt'],
        'in_string_sets': [],
        'not_in_strings': ['spec', 'pectrum', 'structure', 'regno', 'charge', 'over z', 'rsd', 'm/z'],
        'exact_strings': ['m meas.'],
        'type': None,
        'regex': NUMS + r'|' + \
                 LIST_OF_NUMS + r'|' + \
                 NUMS + r'\s*/\s*' + NUMS + r'(\s*Da)?|' + \
                 NUMS + r'\s*-\s*' + NUMS
        },
    
    'parent_moverz_quant': {
        # These are all m/z floating point numbers, but the name is missing the m/z. Talk to Hunter about what to standardize this name to.
        'names': ['parent'],
        'example_values': [],
        'regex_search_strings': [],
        'regex_search_sets': [],
        'not_regex_search_strings': [],
        'in_strings': [],
        'in_string_sets': [],
        'not_in_strings': [],
        'exact_strings': ['parent'],
        'type': 'numeric',
        'regex': None
        },
    'mass_spectrum': {
        # Some of these are the same as just m/z, most look like lists of m/z's or maybe something else.
        # Maybe look to see if the value is just a single m/z and warn about it and change it to m/z.
        # 85:1861.0 86:702.0 87:1148.0   (409.162, 2264.56)(387.1804, 10302.38)
        'regex_search_strings': [],
        'regex_search_sets': [],
        'not_regex_search_strings': [],
        'in_strings': ['spec', 'pectrum'],
        'in_string_sets': [],
        'not_in_strings': ['species', 'composite'],
        'exact_strings': [],
        'type': None,
        'regex': r'(' + NUMS + ':' + NUMS + r'(_|\s+|$))+' + r'|' + \
                 NUMS + r'|' + \
                 LIST_OF_NUMS
        },
    'composite_mass_spectrum': {
        # Unique values when compared to just "spectrum"
        'regex_search_strings': [],
        'regex_search_sets': [],
        'not_regex_search_strings': [],
        'in_strings': [],
        'in_string_sets': [['composite', 'spectrum']],
        'not_in_strings': ['species'],
        'exact_strings': [],
        'type': None,
        'regex': r'(' + PARENTHESIZED_LIST_OF_NUMS + r'\s*)+'
        },
    
    
    'inchi_key': {
        # some values look like 'InChIKey=MTCFGRXMJLQNBG-REOHCLBHSA-N'. standardize to add 'InChIKey='
        # AN004454 has a few rows with 'Internal Standard'. Should this be a special case? I think it should still warn.
        'names': ['inchi key'],
        'regex_search_strings': [],
        'regex_search_sets': [],
        'not_regex_search_strings': [],
        'in_strings': ['inchikey', 'inchi-key', 'inchi_key', 'inchi key'],
        'in_string_sets': [],
        'not_in_strings': [],
        'exact_strings': [],
        'type': None,
        # Some values end with strange characters like '*' or '?' to indicate something about the value, such as it not being confirmed.
        'regex': INCHIKEY + r'(\*?|\?*)' + r'|' +\
                 INCHIKEY_OR_NULL + r'(\s*or\s*|\s*;\s*|\s*_\s*|\s*&\s*)' + INCHIKEY_OR_NULL + r'|' +\
                 r'Sum\s*\(\s*' + INCHIKEY + r'(\s*\+\s*)' + INCHIKEY + r'\s*\)' + r'|' +\
                 LIST_OF_INCHIKEYS_SLASH
        },
    'inchi': {
        # Just a full inchi with InChI= at the start.
        'names': ['inchi', 'base_inchi', 'isotopic_inchi', 'representative_inchi', 'InChI'],
        'regex_search_strings': [],
        'regex_search_sets': [],
        'not_regex_search_strings': [],
        'in_strings': ['inchi'],
        'in_string_sets': [],
        'not_in_strings': ['key'],
        'exact_strings': [],
        'type': None,
        'regex': CHEAP_INCHI + r'|' + BRACKETED_LIST_OF_INCHI
        },
    'smiles': {
        'names': ['SMILES', 'smiles', 'Smiles', 'Canonical SMILES'],
        'example_values': ['C1=CC=C2C(=C1)C(C(=O)N2)CC(=O)O'],
        'regex_search_strings': [],
        'regex_search_sets': [],
        'not_regex_search_strings': [],
        'in_strings': ['smile'],
        'in_string_sets': [],
        'not_in_strings': [],
        'exact_strings': [],
        'type': None,
        'regex': SMILES + r'|' + LIST_OF_SMILES_SEMICOLON
        },
    
    
    'formula': {
        # AN001511 has values like 'C10 H19 N O'. Might want to remove all spaces in this column instead of just stripping.
        # AN002420 has a column 'adducted_isotopic_formula', with values like ['12C12+14N4+16O19+1H32', '12C13+14N5+16O16+1H31+23Na1']
        # Need to be able to handle lists with and without brackets. may need different entry for 'isotopic' formulas with different regex.
        # Formula like 'C9H18NO4+' Ask Hunter about the '+' and adjust regex.  Also 'C5H3H(2)7NO2' and '[13]C3H7[15]NO2' and 'C62H90CoN13O15P-2'
        # The current regex covers a lot, but 'C9H19N3O3･HCl' and 'C62H90CoN13O15P-2' are exceptions. 
        # Note the first one has a wierd character seperating it, it's not a dash like the second one.
        'names': ['formula'],
        'regex_search_strings': [],
        'regex_search_sets': [],
        'not_regex_search_strings': [],
        'in_strings': ['formula'],
        'in_string_sets': [],
        'not_in_strings': [],
        'exact_strings': [],
        'type': 'non-numeric',
        'regex': ISOTOPIC_FORMULA.replace('\d*', '\d*\s*') + r'|' + \
                 CHARGE_FORMULA + r'|' + \
                 GROUP_FORMULA + r'|' + \
                 BRACKETED_LIST_OF_FORMULAS + r'|' + \
                 BRACKETED_LIST_OF_ISOTOPIC_FORMULAS
        },
    
    
    'metabolite': {
        # This is only for validating the required 'Metabolite' column. All other columns with 'metabolite' in the name are handled by other dicts.
        'names': ['metabolite'],
        'regex_search_strings': [],
        'regex_search_sets': [],
        'not_regex_search_strings': [],
        'in_strings': [],
        'in_string_sets': [],
        'not_in_strings': [],
        'exact_strings': ['metabolite'],
        'type': 'non-numeric',
        'regex': None
        },
    'compound': {
        # This is a mixture of IDs, names, and short names. Sort of like metabolite. 
        'names': ['compound', 'Compound Name', 'Compound', 'Compound #', 'Compound.name'],
        'example_values': ['1176.783@0.99987495', '4-Hydroxyethinylestradiol', 'cmp.QI01', 'QI11007', 'X00875', '2HG'],
        'regex_search_strings': [],
        'regex_search_sets': [],
        'not_regex_search_strings': [],
        'in_strings': ['compound', 'compund'],
        'in_string_sets': [],
        'not_in_strings': ['kegg', 'formula', 'pubchem', 'mass', 'rt', 'algo', 'id', 'name'],
        'exact_strings': [],
        'type': 'non-numeric',
        'regex': None
        },
    'name': {
        'names': ['name', ],
        'regex_search_strings': [],
        'regex_search_sets': [],
        'not_regex_search_strings': [],
        'in_strings': ['name'],
        'in_string_sets': [['name', 'refmet']],
        'not_in_strings': ['adduct', 'named', 'internal', 'ion'],
        'exact_strings': [],
        'type': 'non-numeric',
        'regex': None
        },
    'refmet' : {
        'names': [],
        'example_values': [],
        'regex_search_strings': [],
        'regex_search_sets': [],
        'not_regex_search_strings': [],
        'in_strings': ['refmet'],
        'in_string_sets': [],
        'not_in_strings': ['name', 'in'],
        'exact_strings': [],
        'type': None,
        'regex': POSITIVE_INTS
        },
    
    'class': {
        'names': ['class', 'main class', 'sub class', 'super class'],
        'regex_search_strings': [],
        'regex_search_sets': [],
        'not_regex_search_strings': [],
        'in_strings': ['class'],
        'in_string_sets': [],
        'not_in_strings': [],
        'exact_strings': [],
        'type': 'non-numeric',
        'regex': None,
        'inverse_regex': ORGANIC_FORMULA
        },
    'pathway': {
        # Just a string that describes the pathway the metabolite is involved in.
        'names': ['pathway', 'sub pathway', 'super pathway'],
        'example_values': ['amino acids', 'TCA Cycle', 'Endocannabinoid', 'Vitamin A Metabolism'],
        'regex_search_strings': [],
        'regex_search_sets': [],
        'not_regex_search_strings': [],
        'in_strings': ['pathway'],
        'in_string_sets': [],
        'not_in_strings': ['sort'],
        'exact_strings': [],
        'type': 'non-numeric',
        'regex': None
        },
    'pathway_sortorder': {
        # possible column to change to pathway%sortorder
        # Not sure what this is for. Related to pathway somehow, but it's not clear how. Just an integer.
        'names': ['pathway_sortorder'],
        'example_values': ['564'],
        'regex_search_strings': [],
        'regex_search_sets': [],
        'not_regex_search_strings': [],
        'in_strings': [],
        'in_string_sets': [['pathway', 'sort']],
        'not_in_strings': [],
        'exact_strings': [],
        'type': 'integer',
        'regex': None
        },
    
    
    'ion': {
        # There are 2 kinds of columns. Most are numeric m/z columns, some of which will be picked up by the 
        # mass regex. The other kind are columns with values like M+H, -H(-), M+Acid, [M+CH3COO]-. I am not 
        # quite sure what these are about. The difference can only be separated by looking at the values, 
        # not the column names for the most part. Most say 'type' or 'name' or 'species', but some are just labeled 'ion'.
        'regex_search_strings': ['ion', 'ions'],
        'regex_search_sets': [],
        'not_regex_search_strings': [],
        'in_strings': [],
        'in_string_sets': [],
        'not_in_strings': ['adduct', 'm/z', 'mass'],
        'exact_strings': [],
        'type': None,
        'regex': ION + r'|' + \
                 LIST_OF_IONS + r'|' + \
                 LIST_OF_IONS_SPACE + r'|' + \
                 NUMS + r'|' + \
                 LIST_OF_NUMS + r'|' + \
                 NUMS + r'(\s*>\s*|\s*<\s*)' + NUMS + r'|' + \
                 LIST_OF_NUMS_SLASH
        },
    'adduct': {
        # Similar to some of the ion values. Probably more appropriate to say that some ion values were like these values.
        'example_values': ['[M+Cl]-_[M+HCOO]-', '[M+H]+', '[M-H]-', '[M+NH4]+', '[M+H-C5H8O3]+'],
        'regex_search_strings': [],
        'regex_search_sets': [],
        'not_regex_search_strings': [],
        'in_strings': ['adduct'],
        'in_string_sets': [],
        'not_in_strings': ['formula'],
        'exact_strings': [],
        'type': None,
        'regex': ION + r'|' + \
                 LIST_OF_IONS + r'|' + \
                 LIST_OF_IONS_SPACE + r'|' + \
                 LIST_OF_IONS_UNDERSCORE + r'|' + \
                 LIST_OF_IONS_MIXED + r'|' + \
                 LIST_OF_IONS_NO_DELIMITER + r'|' + \
                 BRACKETED_LIST_OF_IONS
        },
    'species': {
        # All values look the same as 'adduct' except for columns 'IS_Species'. Not sure what IS is in 'IS_Species', but the values 
        # are basically abbreviations like LPE(14:0). I'm going to exclude this column for now.
        ## TODO ask Hunter about values like CE(19:0) in the IS_Species column and see what they are and if we can standardize the column name.
        # TODO try to combine the adduct, species, and ion columns that have values like this and make them into 1 name. "ion_species"
        'names': ['Species', 'Ion species', 'Ion Species'],
        'example_values': ['[M+H]+', '[2M+H]+', '[M+Cl]-_[M+CH3COO]-'],
        'regex_search_strings': [],
        'regex_search_sets': [],
        'not_regex_search_strings': [],
        'in_strings': ['species'],
        'in_string_sets': [],
        'not_in_strings': ['is_species', 'ion'],
        'exact_strings': [],
        'type': None,
        # regex needs tested.
        'regex': ION + r'|' + \
                 LIST_OF_IONS + r'|' + \
                 LIST_OF_IONS_UNDERSCORE
        },
    
    
    # IDs
    'pubchem': {
        # other checks, non-negative, integers only
        # Exception values: 'AHLPHDHHMVZTML-BYPYZUCNSA-N', 'HMDB0000092', 'HMDB06029', 'HMDB00746', 'HMDB00619', 'HMDB00653', 'HMDB00563', 'new in this study', '1183?from=summary', '5280880#datasheet=LCSS'
        'names': ['pubchem', 'pubchem id', 'cid'],
        'regex_search_strings': ['cid'],
        'regex_search_sets': [],
        'not_regex_search_strings': [],
        'in_strings': ['pubchem'],
        'in_string_sets': [],
        'not_in_strings': ['formula', 'kegg'],
        'exact_strings': [],
        'type': None,
        # This allows some KEGG values through because it is a common mislabeling.
        'regex': POSITIVE_INTS + r'[&?]?' + r'|' + \
                 LIST_OF_POS_INTS + r'|' + \
                 LIST_OF_POS_INTS_OR + r'|' + \
                 LIST_OF_POS_INTS_SLASH + r'|' + \
                 LIST_OF_POS_INTS_SPACE + r'|' + \
                 LIST_OF_POS_INTS_SEMICOLON + r'|' + \
                 r'Sum \(\d+ \+ \d+\)' + r'|' + \
                 r'CID' + POSITIVE_INTS #+ r'|' + \
                 # r'[CD]\d{5}' + r'|' + \
                 # make_list_regex(r'[CD]\d{5}', r'\|')
        },
    'kegg': {
        # Are IDs that start with other letters besides 'C' okay? Like D01947 is a drug, and ko00260 and CA1373. Probably need to allow '?' characters after IDs.
        # some values look like 'cpd:C00475'. remove the 'cpd:'
        'names': ['kegg', 'kegg id'],
        'regex_search_strings': [],
        'regex_search_sets': [],
        'not_regex_search_strings': [],
        'in_strings': ['kegg'],
        'in_string_sets': [],
        'not_in_strings': ['name'],
        'exact_strings': [],
        'type': None,
        'regex': KEGG + r'|' + \
                 LIST_OF_KEGG + r'|' + \
                 LIST_OF_KEGG_SEMICOLON + r'|' + \
                 LIST_OF_KEGG_SLASH + r'|' + \
                 LIST_OF_KEGG_DOUBLE_SLASH + r'|' + \
                 LIST_OF_KEGG_UNDERSCORE + r'|' + \
                 LIST_OF_KEGG_HYPHEN + r'|' + \
                 LIST_OF_KEGG_MIXED + r'|' + \
                 LIST_OF_KEGG_SPACE + r'|' + \
                 LIST_OF_KEGG_BAR + r'|' + \
                 KEGG + r'-' + FORMULA  + r'|' + \
                 KEGG + r';\d+'
        },
    'hmdb': {
        # Most values are like 'HMDB\d+', but some are just numbers with HMDB dropped. 
        # There are some values with lists of numbers as well as slash separated ones. AN000619
        # Special cas in AN001599, a few values like 'Sum (HMDB11567 + HMDB11537)'
        # There are a few with values like METPA0312 which are MetPA IDs, which are pathway IDs from the MetPA tool, but 
        # it looks to me like that tool is gone now. I'm not sure what to do about these. Talk to Hunter. AN003219
        # There are a few files with columns like 'KEGG or HMDB ID' and 'KEGG_HMDB' which mix IDs. Might need to have a 
        # special case for these kinds and warn about dual use columns as well as expand the values regex to 
        # accept kegg or hmdb values.  AN003618  AN004353
        # exception values 'MXP000013', ', HMDB00510', 'IS', '172.10926', '433.31523', '505.31776', 'EIDNRcv7HeF SpectraBase Compound ID', 'ac107', 'CID57357170', 'ac111'
        'regex_search_strings': [],
        'regex_search_sets': [],
        'not_regex_search_strings': [],
        'in_strings': ['hmdb', 'human metabolome'],
        'in_string_sets': [['hmp', 'id']],
        'not_in_strings': ['class'],
        'exact_strings': [],
        'type': None,
        'regex': HMDB + r'|' + \
                 HMDB_INT + r'|' + \
                 LIST_OF_HMDB + r'|' + \
                 LIST_OF_HMDB_SLASH + r'|' + \
                 LIST_OF_HMDB_INTS + r'|' + \
                 LIST_OF_HMDB_INTS_SLASH + r'|' + \
                 r'Sum \(HMDB\d+ \+ HMDB\d+\)' + r'|' + \
                 LIST_OF_HMDB_AMPERSAND + r'|' + \
                 LIST_OF_HMDB_SEMICOLON + r'|' + \
                 LIST_OF_HMDB_SPACE + r'|' + \
                 r'METPA\d+' + r'|' + \
                 LIST_OF_HMDB_UNDERSCORE
        },
    'lmp': {
        # Lipid Maps ID there are some of these in other_id, should probably detect and warn to make it it's own column.
        # There are a few columns that are "bulk" IDs and I'm not sure exactly how those are different, but they clearly aren't 
        # the same as regular IDs. These seem to have the form 'ST0202-C21H32O2' where it is ID-formula. Not all columns that 
        # have values like that say 'bulk' though. AN002306 has a few values with an asterisk '*' at the end for some reason.
        'example_values': ['LMPK12050102', 'LMST01010204', 'LMGL03010088', 'LMFA01070017', 'ST0202-C21H32O2', 'FA0102-C23H46O2'],
        'regex_search_strings': [],
        'regex_search_sets': [],
        'not_regex_search_strings': [],
        'in_strings': ['lipidmaps', 'lmid'],
        'in_string_sets': [['lmp', 'id'], ['lipid', 'map'], ['lm', 'id']],
        'not_in_strings': [],
        'exact_strings': [],
        'type': None,
        'regex': LIPID_MAPS + r'|' + \
                 LIST_OF_LMP_UNDERSCORE + r'|' + \
                 LIST_OF_LMP + r'|' + \
                 LIST_OF_LMP_SLASH + r'|' + \
                 POSITIVE_INTS
        },
    'chemspider': {
        # Basically the same as PubChem IDs. They are different IDs, but both are just integers with no special form.
        'example_values': [],
        'regex_search_strings': [],
        'regex_search_sets': [],
        'not_regex_search_strings': [],
        'in_strings': ['chemspider'],
        'in_string_sets': [],
        'not_in_strings': [],
        'exact_strings': [],
        'type': None,
        'regex': POSITIVE_INTS + r'[&?]?' + r'|' + \
                 LIST_OF_POS_INTS + r'|' + \
                 LIST_OF_POS_INTS_OR + r'|' + \
                 LIST_OF_POS_INTS_SLASH + r'|' + \
                 LIST_OF_POS_INTS_SPACE + r'|' + \
                 LIST_OF_POS_INTS_SEMICOLON + r'|' + \
                 r'Sum \(\d+ \+ \d+\)' + r'|' + \
                 r'CID' + POSITIVE_INTS + r'|' + \
                 r'CSID' + POSITIVE_INTS + r'|' + \
                 r'[CD]\d{5}' + r'|' + \
                 make_list_regex(r'[CD]\d{5}', r'\|')
        },
    'metlin': {
        # METLIN is a paid database. The IDs seem to just be integers, but there is one dataset like METLIN:21376
        'example_values': ['6888', 'METLIN:21376'],
        'regex_search_strings': [],
        'regex_search_sets': [],
        'not_regex_search_strings': [],
        'in_strings': ['metlin'],
        'in_string_sets': [],
        'not_in_strings': [],
        'exact_strings': [],
        'type': None,
        'regex': POSITIVE_INTS + r'|' + r'METLIN:' + POSITIVE_INTS
        },
    'cas_number': {
        'regex_search_strings': ['cas'],
        'regex_search_sets': [],
        'not_regex_search_strings': [],
        'in_strings': [],
        'in_string_sets': [],
        'not_in_strings': [],
        'exact_strings': [],
        'type': None,
        'regex': CAS + r'|' + LIST_OF_CAS + r'|' + LIST_OF_CAS_SEMICOLON
        },
    'binbase_id': {
        # Must be an integer
        'regex_search_strings': ['bb'],
        'regex_search_sets': [],
        'not_regex_search_strings': [],
        'in_strings': ['binbase'],
        'in_string_sets': [],
        'not_in_strings': [],
        'exact_strings': [],
        'type': None,
        'regex': POSITIVE_INTS
        },
    'chebi_id': {
        # Must be an integer, some values are just integers and some are like CHEBI:1157.
        # The IDs on the ChEBI website are like CHEBI:1157, so I think these should be converted to be like this.
        'regex_search_strings': [],
        'regex_search_sets': [],
        'not_regex_search_strings': [],
        'in_strings': ['chebi'],
        'in_string_sets': [],
        'not_in_strings': [],
        'exact_strings': [],
        'type': None,
        'regex': r'(CHEBI:)?' + POSITIVE_INTS
        },
    'mw_regno': {
        # Identifiers for Metabolomics Workbench structure database.
        'names': ['MW structure'],
        'example_values': ['26', "14991"],
        'regex_search_strings': [],
        'regex_search_sets': [],
        'not_regex_search_strings': [],
        'in_strings': ['regno'],
        'in_string_sets': [['mw', 'structure']],
        'not_in_strings': [],
        'exact_strings': [],
        'type': None,
        'regex': POSITIVE_INTS + r'|' + LIST_OF_POS_INTS
        },
    'mzcloud_id': {
        # ID to the mzCloud database.
        # Integers as you would expect, but some are preceded with 'Reference-' or 'Autoprocessed-'. 
        # This may be meaningful info, but I think it should be removed. Ask Hunter's Opinion.
        'names': ['mzCloud ID', 'mzCloudId'],
        'example_values': ['13', '1744', 'Reference-6554', 'Autoprocessed-8468'],
        'regex_search_strings': [],
        'regex_search_sets': [],
        'not_regex_search_strings': [],
        'in_strings': [],
        'in_string_sets': [['mz', 'cloud', 'id']],
        'not_in_strings': [],
        'exact_strings': [],
        'type': None,
        'regex': r'((Reference|Autoprocessed)-)?' + POSITIVE_INTS
        },
    # Handled under other_id
    # 'Local ID': {
    #     # These are a variety of values from integers to short names to simple row/index numbering.
    #     'regex_search_strings': [],
    #     'regex_search_sets': [],
    #     'not_regex_search_strings': [],
    #     'in_strings': [],
    #     'in_string_sets': [['local', 'id'], ['row', 'id']],
    #     'not_in_strings': [],
    #     'exact_strings': [],
    #     'type': 'ID',
    #     'regex': None
    #     },
    'identifier': {
        # Basically just m/z and retention time mashed together. Not really needed. Should probably print a warning about 
        # making the column 2 separate columns and not having one mashed together.
        # Another variant is mz@rt instead of mz_rt. AN004492 Feature@RT
        # This could just be a subset of local ID.
        'names': ['Identifier', 'identifier'],
        'example_values': ['6.46_600.51_6.45_610.54', '8.45_710.66'],
        'regex_search_strings': [],
        'regex_search_sets': [],
        'not_regex_search_strings': [],
        'in_strings': ['identifier', 'retention time_m/z', 'feature@rt'],
        'in_string_sets': [],
        'not_in_strings': ['pubchem', 'study', 'database'],
        'exact_strings': [],
        'type': None,
        'regex': POS_FLOAT_PAIRS + r'(n|m/z)?' + r'|' + \
                 POS_INT_FLOAT_PAIR + r'|' + \
                 LIST_OF_POS_FLOAT_PAIRS_UNDERSCORE + r'|' + \
                 LIST_OF_POS_FLOAT_PAIRS_NO_SPACE + r'|' + \
                 LIST_OF_POS_FLOAT_PAIRS_MIXED + r'|' + \
                 r'CHEBI:\d+'
        },
    'other_id': {
        # The values vary quite a bit from just integers to traditional DB IDs to a common metabolite name. Can have multiple ID types in one column.
        # For instance, AN000135 has CHEBI:17308 and CID46905267. These seem to come from Pacific Northwest National Laboratory (PNNL_ID).
        # Quite a few labs just set the 'type' as their lab and then the column values are a mix of IDs from databases and other things.
        # A common ID is from Lipid Maps (https://www.lipidmaps.org/databases/lmsd/overview) and their ID is like LMST01020007,
        # LM followed by 2 more capital letters that represent the lipid category and then 8 digits.
        # There is also a LipidMAps Bulk ID (DG (28:0)) and LipidBLAST (DG 27:0) that I don't really understand.
        'example_values': ['CHEBI:17308', 'CID46905267', '747738', '17', 
                           'HMDB00518', 'L-3-Hydroxykynurenine', 'U_Kentucky_isotope',
                           'NST61Q3690', '2Q4710', '12A6030', '16E0950', # AN000202
                           'ADP', 'ADP_13C_1', 'LMST01020007'],
        'regex_search_strings': [],
        'regex_search_sets': [],
        'not_regex_search_strings': ['cas'],
        'in_strings': ['other'],
        'in_string_sets': [['database', 'identifier'], ['chemical', 'id'], ['cmpd', 'id'], 
                           ['database', 'id'], ['database', 'match'], ['local', 'id'], ['row', 'id'],
                           # AN003406 has all 4 of these IDs. Not sure what they go to, but obviously can't all be rename to other_id.
                           # TODO ask Hunter about what to do with multiple other_id columns like this.
                           ['comp', 'id'], ['chem', 'id'], ['chro', 'lib', 'id'], ['lib', 'id']],
        'not_in_strings': ['type', 'pubchem', 'chemspider', 'kegg'],
        'exact_strings': ['id'],
        'type': None,
        'regex': None,
        'inverse_regex': FLOAT + r'|' + SCIENTIFIC_NOTATION
        },
    'other_id_type': {
        # Most types are some string that indicates the ID is a unique one from the lab or university of origin.
        # This could be their own made up system, or a combination of IDs from different sources.
        # Most are also all the same value, but sometimes the type changes to match the value, for instance LipidBLAST and LipidMaps being in the same column.
        # Top used values:
            # LipidMaps                     85
            # UM_Target_ID                  91
            # LipidBLAST                    96
            # LipidMaps Bulk ID            100
        'example_values': ['Osaka_City_University_ID', 'UCDavis_Newman_ID',
                           'HMDB', 'CHEBI', 'LIPIDMAPS', 'LipidBLAST', 'METLIN', 'REFMET', 'BinBase',
                           'UM_Target_ID'],
        'regex_search_strings': [],
        'regex_search_sets': [],
        'not_regex_search_strings': [],
        'in_strings': [],
        'in_string_sets': [['other', 'type'], ['source', 'database']],
        'not_in_strings': [],
        'exact_strings': [],
        'type': 'non-numeric',
        'regex': None
        },
    
        
        
    'retention_time': {
        # Other checks, non-negative, not all zeros (AN004333), has some values with commas instead of decimals. AN003676
        # TODO Need to talk to Hunter about what to do about units, such as (sec) or (min) in column name.
        # Need to look for the columns with values like mz_rt.
        # Some values have a % in them, these are bad and need to be indicated as bad. AN003349
        'names': ['rt', 'retention time'],
        'regex_search_strings': ['rt'],
        'regex_search_sets': [['ret', 'time']],
        'not_regex_search_strings': [],
        'in_strings': ['rtimes', 'r.t.', 'medrt', 'rtsec', 'bestrt', 'compoundrt', 'rtmed'],
        'in_string_sets': [['retention', 'time'], ['rentetion', 'time'], ['retension', 'time']],
        'not_in_strings': ['type', 'error', 'index', 'delta', 'feature', 'm/z'],
        'exact_strings': [],
        'type': None,
        'regex': POSITIVE_FLOATS + r'|' + \
                 r'\d' + r'|' + \
                 POSITIVE_SCIENTIFIC_NOTATION + r'|' + \
                 POSITIVE_FLOAT_RANGE + r'|' + \
                 LIST_OF_POS_FLOATS_UNDERSCORE
        },
    'delta_rt': {
        # MAking delta RT separate because the values can be negative.
        'names': [],
        'regex_search_strings': [],
        'regex_search_sets': [],
        'not_regex_search_strings': [],
        'in_strings': ['deltart'],
        'in_string_sets': [['delta', 'rt']],
        'not_in_strings': ['type', 'error', 'index'],
        'exact_strings': [],
        'type': None,
        'regex': FLOAT + r'|' + r'0'
        },
    'retention_index': {
        # Basically the same as retention time value wise. Mix of integers and floats.
        # Might need a special case to handle names like "Average Rt (min) (LCMS); Retention Index (GCMS)" AN003930
        'names': ['ri', 'retention index'],
        'regex_search_strings': ['ri'],
        'regex_search_sets': [['ret', 'ind'], ['ret', 'index']],
        'not_regex_search_strings': [],
        'in_strings': ['rindex'],
        'in_string_sets': [['retention', 'index'], ['rentetion', 'index'], ['reten', 'index']],
        'not_in_strings': ['type', 'error'],
        'exact_strings': [],
        'type': None,
        'regex': POSITIVE_FLOATS + r'|' + \
                 r'\d' + r'|' + \
                 POSITIVE_SCIENTIFIC_NOTATION + r'|' + \
                 POSITIVE_FLOAT_RANGE + r'|' + \
                 LIST_OF_POS_FLOATS_UNDERSCORE
        },
    'retention_index_type': {
        # Similar to other_id_type. There are a few common index types.
        # Top values:
            # Fiehn               30
            # UFLORIDA_MS_RI      24
            # Mayo_RI             20
            # Time, min           14
            # Nebraska_rt          2
            # HILIC Pos            1
            # RTI_RI               1
            # Jiangnan_ret         1
            # Binbase              1
        'regex_search_strings': [],
        'regex_search_sets': [],
        'not_regex_search_strings': [],
        'in_strings': [],
        'in_string_sets': [['retention', 'index', 'type'], ['ri', 'type']],
        'not_in_strings': ['error'],
        'exact_strings': [],
        'type': 'non-numeric',
        'regex': None
        },
    
    
    # Mostly ignorable columns
    'abbreviation': {
        # Only 2 datasets with this column. One is as you expect, compound name abbreviation. The other looks like a ratio or something.
        'names': ['abbreviation', 'Abbreviation'],
        'example_values': ['aMCA', '12:0'],
        'regex_search_strings': [],
        'regex_search_sets': [],
        'not_regex_search_strings': [],
        'in_strings': ['abbreviation'],
        'in_string_sets': [],
        'not_in_strings': [],
        'exact_strings': [],
        'type': 'non-numeric',
        'regex': None
        },
    'assignment_certainty': {
        # Column giving a numeric value about the confidence in the assigment. Only 3 datasets use it.
        # Numeric values aren't obvious, but might be explained in the mwtabfile or supporting documents.
        'names': [],
        'example_values': ['0', '1', '2'],
        'regex_search_strings': [],
        'regex_search_sets': [],
        'not_regex_search_strings': [],
        'in_strings': [],
        'in_string_sets': [['assignment', 'certainty']],
        'not_in_strings': [],
        'exact_strings': [],
        'type': None,
        'regex': POSITIVE_INTS
        },
    
    
    'comment': {
        # Basically just an arbitrary string value. The couple of columns that have values of 1 and 0 aren't useful and could just be dropped.
        'names': ['Comments', 'Comment', 'comments'],
        'example_values': ['Peak was not detected', '1', '0', 'no MSMS', '2 conformations'],
        'regex_search_strings': [],
        'regex_search_sets': [],
        'not_regex_search_strings': [],
        'in_strings': ['comment'],
        'in_string_sets': [],
        'not_in_strings': [],
        'exact_strings': [],
        'type': None,
        'regex': None,
        'inverse_regex': FLOAT + r'|' + SCIENTIFIC_NOTATION + r'|' + r'\d{2,}'
        },
    'assignment%method': {
        # Column elaborating how the assignment metabolite was assigned. Should normalize the values.
        'names': ['assignment%method'],
        'example_values': ['database', 'indirectly to standard', 'directly to standard', 'direct_to_standard', 'indirect_to_standard'],
        'regex_search_strings': [],
        'regex_search_sets': [],
        'not_regex_search_strings': [],
        'in_strings': [],
        'in_string_sets': [['assignment', 'method']],
        'not_in_strings': [],
        'exact_strings': [],
        'type': 'non-numeric',
        'regex': None
        },
    'isotopologue': {
        # I think these are only in our data sets. Basically just showing how many Carbons, Nitrogens, Hydrogen, etc are the labeled isotope in the compound.
        # The values need some standardization. Need to talk to Hunter about it. I think 'ALL' should probably be removed and the lists should probably lose the brackets.
        'names': ['isotopologue'],
        'example_values': ['13C1', '13C2', 'ALL', "['13C0', '13C0']", 'm0-6', 'm0-m4'],
        'regex_search_strings': [],
        'regex_search_sets': [],
        'not_regex_search_strings': [],
        'in_strings': ['isotopologue'],
        'in_string_sets': [['isotope', 'count']],
        'not_in_strings': ['type'],
        'exact_strings': [],
        'type': 'non-numeric',
        'regex': None
        },
    'isotopologue%type': {
        # Can be gleaned from the isotopologue column, but is a quick way to know exactly which atoms are labeled in the compound.
        'names': ['isotopologue%type'],
        'example_values': ['13C', "['13C+2H', '13C+2H']", "['13C', '13C']", 'C12 PARENT'],
        'regex_search_strings': [],
        'regex_search_sets': [],
        'not_regex_search_strings': [],
        'in_strings': ['isotopologue%type', 'isotope'],
        'in_string_sets': [],
        'not_in_strings': ['count'],
        'exact_strings': [],
        'type': 'non-numeric',
        'regex': None
        },
    'peak_description': {
        # For NMR. More information about the measured peak.
        'names': ['peak_description'],
        'example_values': ["[1H26:C2]HResonance", "['[1H31:13C10]HResonance + [1H31:1H30]J3HH', ' [1H35:13C10]HResonance + [1H35:1H34]J3HH', ' [1H39:13C10]HResonance + [1H39:1H38]J3HH']"],
        'regex_search_strings': [],
        'regex_search_sets': [],
        'not_regex_search_strings': [],
        'in_strings': [],
        'in_string_sets': [['peak', 'description']],
        'not_in_strings': [],
        'exact_strings': [],
        'type': 'non-numeric',
        'regex': None
        },
    'peak_pattern': {
        # For NMR. Short name describing the peak pattern, such as a single peak or 2 or 3 peaks.
        # TODO Talk to Hunter about some standardization of values. Particularly removing the parantheses and asking about single chars vs full words.
        'names': ['peak_pattern', 'Peak patterns'],
        'example_values': ['multiplets', 'singlet', 'ddd', 'Singlet', 'doublets', 'triplet',
                           'q', 'dd', '(t)', '(d)', '(m)', '(s)', '(dd)', 'Triplet', 'Quatet',
                           'multiplet', 'quatet', 'Doublets'],
        'regex_search_strings': [],
        'regex_search_sets': [],
        'not_regex_search_strings': [],
        'in_strings': [],
        'in_string_sets': [['peak', 'pattern']],
        'not_in_strings': [],
        'exact_strings': [],
        'type': 'non-numeric',
        'regex': None
        },
    'transient_peak': {
        # For NMR. I can't quite remember what this is about, I think it is numbering the peaks used for 1 compound.
        'names': ['transient_peak'],
        'example_values': ['1', '2', '3', '4', '5', '6'],
        'regex_search_strings': [],
        'regex_search_sets': [],
        'not_regex_search_strings': [],
        'in_strings': [],
        'in_string_sets': [['transient', 'peak']],
        'not_in_strings': ['type'],
        'exact_strings': [],
        'type': None,
        'regex': POSITIVE_INTS
        },
    'transient_peak%type': {
        # For NMR. There is only 1 value. Might want to talk to Hunter about what the transient_peak stuff is again.
        'names': ['transient_peak%type'],
        'example_values': ['integer'],
        'regex_search_strings': [],
        'regex_search_sets': [],
        'not_regex_search_strings': [],
        'in_strings': [],
        'in_string_sets': [['transient', 'peak', 'type']],
        'not_in_strings': [],
        'exact_strings': [],
        'type': 'non-numeric',
        'regex': None
        },
    'FISh Coverage': {
        # "FISh" coverage mass spectrometry refers to a technique in proteomics where "FISh" stands for 
        # "Full-length Isotope-labeled Stable Heavy peptide" and is used to measure the extent of protein 
        # sequence coverage achieved in a mass spectrometry analysis by comparing the detected peptide 
        # fragments to a reference protein sequence
        # This seems to be specific to a few data sets. and is only numeric.
        'names': ['FISh.Coverage'],
        'example_values': ['10', "38.89"],
        'regex_search_strings': [],
        'regex_search_sets': [],
        'not_regex_search_strings': [],
        'in_strings': [],
        'in_string_sets': [['fish', 'coverage']],
        'not_in_strings': [],
        'exact_strings': [],
        'type': None,
        'regex': POSITIVE_NUMS
        },
    'MSI Category': {
        # MSI seems to stand for Mass Spec Imaging, but not sure what "category" is. The only value is 1 across all datasets.
        'names': ['MSIcategory'],
        'example_values': ['1'],
        'regex_search_strings': [],
        'regex_search_sets': [],
        'not_regex_search_strings': [],
        'in_strings': ['msicategory'],
        'in_string_sets': [],
        'not_in_strings': [],
        'exact_strings': [],
        'type': None,
        'regex': r'1'
        },
    'annotations': {
        # String that has more information like IDs all in 1 string. Sometiems just a string that is a name.
        # Ask Hunter about warning on this type of column to create actual columns for the data. These should have their own columns.
        # Also ask about creating the columns in repair command.
        'names': ['annotation', 'annotations', 'detailed annotation'],
        'example_values': ['[ C6 H3 O5, mfg=75.14, overall=75.14 ]', 
                           'Creatine [ C4 H9 N3 O2, mfg=99.07, db=86.69, tgt=, overall=92.88, CAS ID=57-00-1, KEGG ID=C00300, HMP ID=HMDB00064, METLIN ID=7 ]',
                           'Cer(d18:1/14:0)',
                           'Analyte List - Metabolite'],
        'regex_search_strings': [],
        'regex_search_sets': [],
        'not_regex_search_strings': ['id'],
        'in_strings': ['annotation'],
        'in_string_sets': [],
        'not_in_strings': ['source', 'approach', 'confidence', 'level'],
        'exact_strings': [],
        'type': 'non-numeric',
        'regex': None
        },
    'ISTD': {
        # string of a name. I think this might be internal standard.
        'names': ['ISTD'],
        'example_values': ['1_CE (22:1)', '1_Cer (d18:1/17:0)', '1_DG (12:0/12:0/0:0)'],
        'regex_search_strings': [],
        'regex_search_sets': [],
        'not_regex_search_strings': [],
        'in_strings': ['internal'],
        'in_string_sets': [],
        'not_in_strings': [],
        'exact_strings': ['istd'],
        'type': 'non-numeric',
        'regex': None
        },
    
    
    
    # POS and NEG
    'platform': {
        # Mostly trying to indicate positive or negative mode, LC or GC. I think these values could be standardized, talk to Hunter.
        # similar to polarity, esi, ionization, should try to combine and or separate these.
        # TODO
        'names': ['platform', 'PLATFORM'],
        'example_values': ['LC/MS Neg', 'LC/MS Polar', 'LC/MS Pos Late',
                           'LC/MS Pos Early', 'Pos Early',
                           'Polar', 'Neg', 'LC/MS Pos',
                           'GC/MS', 'Pos Late'],
        'regex_search_strings': [],
        'regex_search_sets': [],
        'not_regex_search_strings': [],
        'in_strings': ['platform'],
        'in_string_sets': [],
        'not_in_strings': [],
        'exact_strings': [],
        'type': 'non-numeric',
        'regex': None
        },
    'MS Method': {
        # String describing the mass spec method, but some values seem more like a class of compounds than a method to me.
        'names': ['Method', 'MS_method'],
        'example_values': ['HILIC_QE_Pos', 'HILIC_QE_Neg', 'HILIC_TQS_Pos', 'HILIC_TQS_Neg', 'HILIC-pos',
                          'HILIC-neg', 'C18-neg', 'C8-pos', 'Primary Metabolomics', 'Biogenic Amines', 'Lipidomics', 'RP_TQS_Pos'],
        'regex_search_strings': [],
        'regex_search_sets': [],
        'not_regex_search_strings': [],
        'in_strings': ['method'],
        'in_string_sets': [],
        'not_in_strings': ['assignment'],
        'exact_strings': [],
        'type': 'non-numeric',
        'regex': None
        },
    'polarity': {
        # Basically just positive or negative, but there are variations on the values, such as pos, positive, and +.
        # There are a few datasets that look like they have this mislabeled and it should be adduct or ion or something like that. Should probably discuss with Hunter.
        # We can likely standardize values.
        'names': ['Polarity', 'polarity', 'MS polarity', 'POLARITY'],
        'example_values': ['Negative', 'Positive', 'pos', 'neg', 'positive', 'negative', '+', '-1', '1', 'POSITIVE', 'NEGATIVE'],
        'regex_search_strings': [],
        'regex_search_sets': [],
        'not_regex_search_strings': [],
        'in_strings': ['polarity'],
        'in_string_sets': [],
        'not_in_strings': [],
        'exact_strings': [],
        'type': None,
        'regex': r'(?i)(neg|pos|1|-1|\+|positive|negative)' + r'|\[M\+H\]\+|\[M-H\]-|5MM\+|5MM-'
        },
    'ESI mode': {
        # Very similar to 'polarity'. Can probably standardize values, and should look into whether this and polarity are just synonyms.
        # May also be the same as positive or negative MS mode, and should be unnecessary since this info should be in the MS section.
        # Although for ESI columns there are mixed positive and negative, so this might not be the case. Or is a bad data set.
        'names': ['ESI mode', 'ESI'],
        'example_values': ['ESI (+)', '(-)', 'ESI (-)', '(-) ES)', '-1', '(+) ESI', '(-) ESI', 'neg', 'pos'],
        'regex_search_strings': [],
        'regex_search_sets': [],
        'not_regex_search_strings': [],
        'in_strings': ['esi'],
        'in_string_sets': [],
        'not_in_strings': [],
        'exact_strings': [],
        'type': None,
        'regex': r'(neg|pos|1|-1|(ESI )?\(\+\)( ESI)?|(ESI )?\(-\)( ESI| ES\))?|positive|negative)'
        },
    'Ionization mode': {
        # Very similar to 'polarity' and 'ESI mode'. One column even has values like 'ES-'
        'names': ['Ionization mode', 'used ionization mode'],
        'example_values': ['+', 'ES-', 'ES+', 'Positive', 'Negative', 'positive', 'negative'],
        'regex_search_strings': [],
        'regex_search_sets': [['pos', 'neg']],
        'not_regex_search_strings': [],
        'in_strings': ['ionization', 'ionisation'],
        'in_string_sets': [],
        'not_in_strings': ['confirmed'],
        # 5 out of 6 columns named mode have values that are Positive or Negative, the other one is something else.
        # 2 out of 4 columns named MS mode have values like pos and neg. The other is TOF.
        'exact_strings': ['mode', 'ms mode'],
        'type': None,
        'regex': r'(?i)(neg|pos|1|-1|(ES)?\+|(ES)?-|positive|negative|TOF|Splitless|Split30)'
        },
    
    
    
    'frequency': {
        # Just an integer. I think it's like how many scans that m/z or compound appeared in.
        'names': ['frequency'],
        'example_values': ['32', '0', '110'],
        'regex_search_strings': [],
        'regex_search_sets': [],
        'not_regex_search_strings': [],
        'in_strings': ['frequency'],
        'in_string_sets': [],
        'not_in_strings': [],
        'exact_strings': [],
        'type': None,
        'regex': POSITIVE_INTS
        },
    }


for key, attributes in column_matching_attributes.items():
    # The extra paranthese around the regular expressions are needed to use the PyArrow backend in pandas.
    # As of 3-12-2025 this is required, but there is an open issue (https://github.com/pandas-dev/pandas/issues/61072) that might get resolved later.
    column_matching_attributes[key]['values_regex'] = r'(' + attributes['regex'] + r')' if attributes['regex'] else None
    column_matching_attributes[key]['values_inverse_regex'] = r'(' + attributes['inverse_regex'] + r')' if 'inverse_regex' in attributes else None







