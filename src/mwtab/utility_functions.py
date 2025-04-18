# -*- coding: utf-8 -*-
"""
Various functions used in more than one module.
"""
import re
import itertools

import pandas


def _create_numeric_df(df):
    """
    """
    numeric_df = df.iloc[:, 1:]
    numeric_df = numeric_df.apply(pandas.to_numeric, errors='coerce')
    # 0 should be the same as nan for this, AN000027.
    numeric_df[numeric_df.isna()] = 0
    return numeric_df


FAMILY_REGEX = r'(_(\d+)|-ALL|-13C(\d+)|-15N(\d+)|-[13C_15N](\d+))'
def find_metabolite_families(metabolites):
    """
    """
    roots = []
    for met in metabolites:
        if re_match := re.fullmatch(r'(.*)' + FAMILY_REGEX, met):
            roots.append(re_match.group(1))
        else:
            roots.append(met)
    roots = list(set(roots))
    root_to_mets = {root:[root] for root in roots}
    met_to_root = {root:root for root in roots}
    for met in metabolites:
        if met in roots:
            continue
        for root in roots:
            if re.fullmatch(re.escape(root) + FAMILY_REGEX, met):
                met_to_root[met] = root
                root_to_mets[root].append(met)
    
    return met_to_root, root_to_mets

def find_family_sequences(root_to_mets):
    """
    """
    root_to_sequence = {}
    for root, mets in root_to_mets.items():
        if len(mets) == 1:
            continue
        sequence = []
        for met in mets:
            if re_match := re.fullmatch(r'.*' + FAMILY_REGEX, met):
                sequence_num = [group for group in re_match.groups()[1:] if group and re.fullmatch(r'\d+', group)]
                if sequence_num:
                    sequence.append(int(sequence_num[0]))
        if sequence:
            root_to_sequence[root] = sorted(sequence)
    return root_to_sequence


def _determine_row_subsets(df):
    """
    """
    row_sets = [set([value for value in row.values if not pandas.isna(value)]) for index, row in df.iterrows()]
    is_a_subset = [False]*df.shape[0]
    for perm in itertools.permutations(range(df.shape[0]), r=2):
        is_a_subset[perm[0]] |= row_sets[perm[0]].issubset(row_sets[perm[1]])
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








