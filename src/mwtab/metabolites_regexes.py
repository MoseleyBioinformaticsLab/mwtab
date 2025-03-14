# -*- coding: utf-8 -*-
"""
Regular expressions to match values in METABOLITES blocks.
"""



def make_list_regex(element_regex, delimiter, quoted_elements=False, empty_string=False):
    """ Creates a regular expression that will match a list of element_regex delimited by delimiter.
    
    Note that delimiter can be a regular expression like (,|;) to match 2 different types of delimiters. 
    If quoated_elements is True, then allow element_regex to surrounded by single or double quotes. 
    If empty_string is True, then the list regex will match a single element_regex and the empty string.
    """
    if quoted_elements:
        # The simplest regex to add quotes around the element, but allows elements like 'element" through.
        # element_regex = r'(\'|"|)' + element_regex + r'(\'|"|)'
        # Forces an element to be surrounded by the same thing on both sides, but still allows for mixed elements.
        # element_regex = r"('" + element_regex + r"'" + r'|' + r'"' + element_regex + r'"' + r'|' + element_regex + r')'
        
        # All elements must have the same type of quotation marks.
        return r'(' + make_list_regex(element_regex, delimiter, empty_string=empty_string) + r'|' + \
               make_list_regex(f"'{element_regex}'", delimiter, empty_string=empty_string) + r'|' + \
               make_list_regex(f'"{element_regex}"', delimiter, empty_string=empty_string) + r')'
    
    repetition_symbol = '*' if empty_string else '+'
    
    return r'((' + element_regex + r'\s*' + delimiter + r'\s*)' + repetition_symbol + r'(' + element_regex + r'\s*|\s*))'



# # Testing make_list_regex
# pass_values = ['a,', 'a,a', 'a,a,', 'a ,', 'a , a', 'a , a ,']
# fail_values = ['a', '']

# # pass should pass and fail should fail
# regex = make_list_regex('a', ',')
# # pass should pass, but fail should also pass now.
# regex = make_list_regex('a', ',', False, True)

# for value in pass_values:
#     if not re.fullmatch(regex, value):
#         print(value)

# for i, value in enumerate(fail_values):
#     if re.fullmatch(regex, value):
#         print(i, value)


# pass_values = ['["a",]', '["a","a"]', '["a","a",]', '["a" ,]', '["a" , "a"]', '["a" , "a" ,]', '["a"]', '[]']
# fail_values = ['']

# regex = r'\[' + make_list_regex('a', ',', True, True) + r'\]'



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


