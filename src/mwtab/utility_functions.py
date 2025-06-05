# -*- coding: utf-8 -*-
"""
Various functions used in more than one module.
"""
import re
import itertools

import pandas
import fuzzywuzzy.fuzz


def _create_numeric_df(df):
    """
    """
    numeric_df = df.drop(['Metabolite'], axis=1)
    numeric_df = numeric_df.apply(pandas.to_numeric, errors='coerce')
    # 0 should be the same as nan for this, AN000027.
    numeric_df[numeric_df.isna()] = 0
    return numeric_df


# FAMILY_REGEX = r'((_(\d+))+|-ALL|-13C(\d+)|-15N(\d+)|-[13C_15N](\d+))'
def find_metabolite_families(metabolites):
    """Look for familes like '6-Keto Prostaglandin F1', '6-Keto Prostaglandin F1_1', '6-Keto Prostaglandin F1_2' in metabolites.
    
    Families can also include sets like ['root_name', 'root_name-ALL', 'root_name-15N0']. Go over the 
    list of metabolites and identify roots and metabolites. Group them together into 2 dictionaries, 
    met_to_root which is a 1-1 mapping of metabolites to their root, which could be themselves. root_to_mets, 
    which is a maapping of roots to a list of all metabolites that share that root. The list includes 
    the root itself. 
    
    For families that have names ending in (_\d+) also create a dictionary, root_to_sequence, 
    which maps roots to a sorted list of integers. For example, the family 
    ['6-Keto Prostaglandin F1', '6-Keto Prostaglandin F1_1', '6-Keto Prostaglandin F1_2'] would result in a 
    root_to_sequence of {'6-Keto Prostaglandin F1': [1, 2]}. Note that root_to_sequence identifies 
    roots slightly differently. For example, the family ['LPC', 'LPC_16_0', 'LPC_17_6', 'LPC_18_2', 'LPC_18_1', 'LPC_18_0', 'LPC_19_6'] 
    will generate a root for 'LPC_18' and nothing else, but root_to_mets will have a root for 'LPC' and nothing else.
    
    The slight differences in roots is due to how the dictionaries are used later, but they are all 
    generated in this function because of a certain amount of shared code and the slow nature of regular 
    expressions needing optimization.
    
    Args:
        metabolites: list of names of metabolites
    
    Returns:
        3 dictionaries, met_to_root, root_to_mets, root_to_sequence described in the function description.
    """
    family_regex = r'^(.*?)' + r'((_(\d+))+|-ALL|-13C(\d+)|-15N(\d+)|-[13C_15N](\d+))' + r'$'
    met_series = pandas.Series(metabolites, dtype='string[pyarrow]')
    tail_groups = met_series.str.extract(family_regex)
    roots = tail_groups.iloc[:, 0]
    na_roots_mask = roots.isna()
    roots.loc[na_roots_mask] = met_series[na_roots_mask]
    roots_to_match = roots[~na_roots_mask]
    root_groups = roots_to_match.groupby(roots_to_match)
    root_to_mets = {met:[met] for met in roots[na_roots_mask]}
    met_to_root = {met:met for met in metabolites}
    for name, root_group in root_groups:
        root = root_group.iloc[0]
        mets = [root] + list(met_series.loc[root_group.index])
        root_to_mets[root] = mets
        met_to_root.update({met:root for met in mets})
    
    # Determine the sequence for roots where that's relevant.
    root_to_sequence = {}
    sequence_roots_mask = tail_groups.iloc[:, 1].str.fullmatch(r'(_(\d+))+', na=False)
    if sequence_roots_mask.any():
        # multi_num_roots_mask = tail_groups.iloc[:, 1].str.fullmatch(r'(_(\d+)){2,}', na=False)
        sequence_tail_groups = tail_groups.loc[sequence_roots_mask, :]
        sequence_split_df = sequence_tail_groups.iloc[:, 1].str.split('_', regex=False)
        sequence_str_to_concat = sequence_split_df.apply(lambda x: '_'.join(x[0:-1]))
        sequence_roots = sequence_tail_groups.iloc[:, 0] + sequence_str_to_concat
        sequence_roots = sequence_roots.loc[sequence_roots.duplicated(keep=False)]
        # sequence_nums = sequence_tail_groups.loc[sequence_roots.index, 3]
        sequence_groups = sequence_roots.groupby(sequence_roots)
        for name, sequence_group in sequence_groups:
            root = sequence_group.iloc[0]
            sequence = list(sequence_tail_groups.loc[sequence_group.index, 3].astype(int))
            root_to_sequence[root] = sorted(sequence)
    
    return met_to_root, root_to_mets, root_to_sequence

# def find_family_sequences(root_to_mets):
#     """
#     """
#     root_to_sequence = {}
#     for root, mets in root_to_mets.items():
#         if len(mets) == 1:
#             continue
#         sequence = []
#         for met in mets:
#             if re_match := re.fullmatch(r'.*' + FAMILY_REGEX, met):
#                 sequence_num = [group for group in re_match.groups()[1:] if group and re.fullmatch(r'\d+', group)]
#                 if sequence_num:
#                     sequence.append(int(sequence_num[0]))
#         if sequence:
#             root_to_sequence[root] = sorted(sequence)
#     return root_to_sequence

def _determine_row_subsets(df):
    """
    """
    row_sets = [set([(col, value) for col, value in row.items() if not pandas.isna(value)]) for index, row in df.iterrows()]
    is_a_subset = [False]*df.shape[0]
    for perm in itertools.permutations(range(df.shape[0]), r=2):
        is_a_subset[perm[0]] |= (row_sets[perm[0]].issubset(row_sets[perm[1]]) and row_sets[perm[0]] != row_sets[perm[1]])
    return is_a_subset

def drop_duplicate_subsets(df):
    """
    Note that duplicate rows will drop both rows, so make sure to run df.drop_duplicates() 
    before this if you don't want that.
    """
    metabolites = df.loc[:, ['Metabolite']]
    duplicates_bool = metabolites.duplicated(keep=False)
    if duplicates_bool.empty:
        return df
    duplicates = df[duplicates_bool]
    groups = duplicates.groupby(['Metabolite'], dropna=False)
    indexes_to_drop = []
    for name, group in groups:
        is_a_subset = _determine_row_subsets(group)
        # row_sets = [set([value for value in row.values if not pandas.isna(value)]) for index, row in group.iterrows()]
        # is_a_subset = [False]*group.shape[0]
        # for perm in itertools.permutations(range(group.shape[0]), r=2):
        #     is_a_subset[perm[0]] |= row_sets[perm[0]].issubset(row_sets[perm[1]])
        indexes_to_drop += list(group.index[[i for i, is_subset in enumerate(is_a_subset) if is_subset]])
    return df.loc[~df.index.isin(indexes_to_drop), :]


def _determine_permutation_products_of_pairs(character_set_pairs: list[list[list[str]]]):
    """
    """
    character_set_products = []
    for character_set in character_set_pairs:
        perms = [[perm for perm in itertools.permutations(character_pair, r=2)] for character_pair in character_set]
        character_set_products.append([prod for prod in itertools.product(*perms)])
    return character_set_products

CHARACTER_SETS_FOR_REPLACEMENT = [[['^', "'"], ['_', ',']], [['+', '_']], [[' + ', '_']]]
character_set_products = _determine_permutation_products_of_pairs(CHARACTER_SETS_FOR_REPLACEMENT)
REPLACEMENT_SETS = [item for sublist in character_set_products for item in sublist]

def _chain_replace(list_of_pairs, string):
    """
    """
    for pair in list_of_pairs:
        string = string.replace(pair[0], pair[1])
    return string

def _compute_fuzz_ratio(string1, string2, inclusion_ratio=90, replacement_sets = REPLACEMENT_SETS):
    """
    """
    lowered_string1 = string1.lower()
    lowered_string2 = string2.lower()
    if len(string1) == len(string2):
        replacements = [_chain_replace(replacement_set, lowered_string1) for replacement_set in REPLACEMENT_SETS]
        if any([replacement == lowered_string2 for replacement in replacements]):
            return 100
        else:
            for i in range(len(string1)):
                replace_char1 = None
                if lowered_string1[i] != lowered_string2[i]:
                    replace_char1 = lowered_string1[i]
                    replace_char2 = lowered_string2[i]
                    break
            # Don't replace characters that are numbers or letters, only symbols.
            # Otherwise 'VITAMIN D3_1' and 'VITAMIN D3_2' get marked as the same. AN000063
            if replace_char1 and \
               not replace_char1.isalnum() and \
               not replace_char2.isalnum() and \
               lowered_string1.replace(replace_char1, replace_char2) == lowered_string2:
                return 100
    
    ratio = fuzzywuzzy.fuzz.ratio(lowered_string1, lowered_string2)
    if ratio >= inclusion_ratio:
        match1 = re.fullmatch(r'(.*)_(\d+)', lowered_string1)
        match2 = re.fullmatch(r'(.*)_(\d+)', lowered_string2)
        if match1 and match2 and match1.group(2) != match2.group(2):
            return 0
        # Don't want 'VITAMIN D3' fuzzy matching to 'VITAMIN D3_1'.
        elif match1:
            replacements = [_chain_replace(replacement_set, lowered_string2) for replacement_set in REPLACEMENT_SETS]
            if lowered_string2 in lowered_string1 or any([replacement in lowered_string1 for replacement in replacements]):
                return 0
        elif match2:
            replacements = [_chain_replace(replacement_set, lowered_string1) for replacement_set in REPLACEMENT_SETS]
            if lowered_string2 in lowered_string1 or any([replacement in lowered_string2 for replacement in replacements]):
                return 0
        return ratio

def compute_fuzz_ratios(list1, list2, inclusion_ratio=90):
    """
    """
    fuzz_ratios = {}
    for element1 in list1:
        temp_dict = {}
        for element2 in list2:
            ratio = _compute_fuzz_ratio(element1, element2, inclusion_ratio)
            if ratio:
                temp_dict[element2] = ratio
        
        if temp_dict:
            fuzz_ratios[element1] = pandas.Series(temp_dict).sort_values()
    return fuzz_ratios





