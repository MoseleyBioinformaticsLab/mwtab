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


FAMILY_REGEX = r'((_(\d+))+|-ALL|-13C(\d+)|-15N(\d+)|-[13C_15N](\d+))'
def find_metabolite_families(metabolites):
    """
    """
    family_regex = r'^(.*?)' + FAMILY_REGEX + r'$'
    met_series = pandas.Series(metabolites, dtype='string[pyarrow]')
    roots = met_series.str.extract(family_regex).iloc[:, 0]
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
    return met_to_root, root_to_mets

# TODO modify this to handle groups like ['LPC', 'LPC_16_0', 'LPC_17_6', 'LPC_18_2', 'LPC_18_1', 'LPC_18_0', 'LPC_19_6'] from AN004198.
# Currently, only the last integer will be used to sequence instead of the whole thing. Need to look for matches 
# like _18_2 and transform to 182 for sequencing.
def find_family_sequences(root_to_mets):
    """
    """
    root_to_sequence = {}
    for root, mets in root_to_mets.items():
        if len(mets) == 1:
            continue
        sequence = []
        for met in mets:
            if re_match := re.fullmatch(r'.*?' + FAMILY_REGEX, met):
                sequence_num = [group for group in re_match.groups()[1:] if group and re.fullmatch(r'\d+', group)]
                if sequence_num:
                    sequence.append(int(sequence_num[0]))
        if sequence:
            root_to_sequence[root] = sorted(sequence)
    return root_to_sequence

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


def _compute_fuzz_ratio(string1, string2, inclusion_ratio=90):
    """
    """
    lowered_string1 = string1.lower()
    lowered_string2 = string2.lower()
    if len(string1) == len(string2):
        replacement1 = lowered_string1.replace('^', "'")
        replacement1 = replacement1.replace('_', ',')
        replacement2 = lowered_string1.replace("'", "^")
        replacement2 = replacement2.replace(',', '_')
        replacement3 = lowered_string1.replace("'", "^")
        replacement3 = replacement3.replace('_', ',')
        replacement4 = lowered_string1.replace("^", "'")
        replacement4 = replacement4.replace(',', '_')
        replacement5 = lowered_string1.replace('+', '_')
        replacement6 = lowered_string1.replace('_', '+')
        if replacement1 == lowered_string2 or \
           replacement2 == lowered_string2 or \
           replacement3 == lowered_string2 or \
           replacement4 == lowered_string2 or \
           replacement5 == lowered_string2 or \
           replacement6 == lowered_string2:
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
            replacement1 = lowered_string2.replace('^', "'")
            replacement1 = replacement1.replace('_', ',')
            replacement2 = lowered_string2.replace("'", "^")
            replacement2 = replacement2.replace(',', '_')
            replacement3 = lowered_string2.replace("'", "^")
            replacement3 = replacement3.replace('_', ',')
            replacement4 = lowered_string2.replace("^", "'")
            replacement4 = replacement4.replace(',', '_')
            replacement5 = lowered_string2.replace('+', '_')
            replacement6 = lowered_string2.replace('_', '+')
            if lowered_string2 in lowered_string1 or \
               replacement1 in lowered_string1 or \
               replacement2 in lowered_string1 or \
               replacement3 in lowered_string1 or \
               replacement4 in lowered_string1 or \
               replacement5 in lowered_string1 or \
               replacement6 in lowered_string1:
                return 0
        elif match2:
            replacement1 = lowered_string1.replace('^', "'")
            replacement1 = replacement1.replace('_', ',')
            replacement2 = lowered_string1.replace("'", "^")
            replacement2 = replacement2.replace(',', '_')
            replacement3 = lowered_string1.replace("'", "^")
            replacement3 = replacement3.replace('_', ',')
            replacement4 = lowered_string1.replace("^", "'")
            replacement4 = replacement4.replace(',', '_')
            replacement5 = lowered_string1.replace('+', '_')
            replacement6 = lowered_string1.replace('_', '+')
            if lowered_string1 in lowered_string2 or \
               replacement1 in lowered_string2 or \
               replacement2 in lowered_string2 or \
               replacement3 in lowered_string2 or \
               replacement4 in lowered_string2 or \
               replacement5 in lowered_string2 or \
               replacement6 in lowered_string2:
                return 0
        return ratio

def compute_fuzz_ratios(list1, list2):
    """
    """
    fuzz_ratios = {}
    for element1 in list1:
        temp_dict = {}
        for element2 in list2:
            ratio = _compute_fuzz_ratio(element1, element2)
            if ratio:
                temp_dict[element2] = ratio
        
        if temp_dict:
            fuzz_ratios[element1] = pandas.Series(temp_dict).sort_values()
    return fuzz_ratios





