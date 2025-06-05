# -*- coding: utf-8 -*-
"""
Augment the DATA and METABOLITES tables of an mwtab file with extra columns.
"""
import re
import json
import logging
from contextlib import redirect_stdout, redirect_stderr
import datetime

import pandas

from . import utility_functions
from . import repair_metabolites_matching

_create_numeric_df = utility_functions._create_numeric_df
find_metabolite_families = utility_functions.find_metabolite_families
_determine_row_subsets = utility_functions._determine_row_subsets
_compute_fuzz_ratio = utility_functions._compute_fuzz_ratio

column_finders = repair_metabolites_matching.column_finders


def get_duplicate_rows(df, section_name, numeric=True):
    """
    """    
    # Convert every column except the metabolite column to numbers.
    if numeric:
        numeric_df = _create_numeric_df(df)
    else:
        numeric_df = df.drop(['Metabolite'], axis=1)
    
    duplicates_bool = numeric_df.duplicated(keep=False)
    if duplicates_bool.empty:
        return pandas.DataFrame()
    duplicates = numeric_df[duplicates_bool]
    # Remove rows that are mono value. For example, if a row has all nan's that should not be considered as a duplicate.
    duplicates = duplicates[~duplicates.nunique(axis = 1, dropna=False).eq(1)]
    
    return duplicates


def _count_from_match(re_match):
    """
    """
    if re_match:
        if re_match.group(2):
            group_num = 3
        elif re_match.group(5):
            group_num = 6
        count = 0 if not re_match.group(group_num) else int(re_match.group(group_num))
    else:
        count = None
    
    return count

def get_isotope_counts(iso_str):
    """
    """
    iso_str = iso_str.lower()
    carbon_match = re.fullmatch(r'.*(?:[^\d]|^|2h\d*|15n\d*)((13c(?:_|\-)?(\d*))(15n|2h)|(13c(?:_|\-)?(\d*))).*', iso_str)
    nitrogen_match = re.fullmatch(r'.*(?:[^\d]|^|2h\d*|13c\d*)((15n(?:_|\-)?(\d*))(13c|2h)|(15n(?:_|\-)?(\d*))).*', iso_str)
    hydrogen_match = re.fullmatch(r'.*(?:[^\d]|^|13c\d*|15n\d*)((2h(?:_|\-)?(\d*))(15n|13c)|(2h(?:_|\-)?(\d*))).*', iso_str)
    
    carbon_count = _count_from_match(carbon_match)
    nitrogen_count = _count_from_match(nitrogen_match)
    hydrogen_count = _count_from_match(hydrogen_match)
    
    # Only change None to 0 if at least 1 count is found.
    if carbon_count is not None or nitrogen_count is not None or hydrogen_count is not None:
        carbon_count = carbon_count if carbon_count is not None else 0
        nitrogen_count = nitrogen_count if nitrogen_count is not None else 0
        hydrogen_count = hydrogen_count if hydrogen_count is not None else 0
    
    return carbon_count, nitrogen_count, hydrogen_count
    
def remove_isotope_substring(compound):
    """
    """
    original_compound = compound
    compound = re.sub(r'([^\d]|^|2h\d*|15n\d*)(?:(?:13c(?:_|\-)?(?:\d*))(15n|2h)|(?:13c(?:_|\-)?(?:\d*)))', _handle_replacement, compound)
    compound = re.sub(r'([^\d]|^|2h\d*|13c\d*)(?:(?:15n(?:_|\-)?(?:\d*))(13c|2h)|(?:15n(?:_|\-)?(?:\d*)))', _handle_replacement, compound)
    compound = re.sub(r'([^\d]|^|13c\d*|15n\d*)(?:(?:2h(?:_|\-)?(?:\d*))(15n|13c)|(?:2h(?:_|\-)?(?:\d*)))', _handle_replacement, compound)
    if compound != original_compound and original_compound[-1] != '-':
        compound = compound.strip('_- ')
    return compound

def _handle_replacement(match_object):
    """
    """
    if match_object.group(2):
        return match_object.group(1) + match_object.group(2)
    return match_object.group(1)


def fuzzy_match(working_row, df, max_len_diff):
    """
    """
    working_df = df[~df.index.isin([working_row.name])]
    working_len = len(working_row['lowered_name'])
    working_row_has_lipid_nums = not pandas.isna(working_row['lipid_nums'])
    ratios = {}
    for index, row in working_df.iterrows():
        row_len = len(row['lowered_name'])
        len_diff = abs( working_len - row_len)
        row_has_lipid_nums = not pandas.isna(row['lipid_nums'])
        
        if (working_row_has_lipid_nums and not row_has_lipid_nums) or (row_has_lipid_nums and not working_row_has_lipid_nums):
            ratio = 0
        elif working_row['lowered_name'] in row['lowered_name'] and len_diff < max_len_diff and len_diff < min(row_len, working_len):
            ratio = 100
        else:
            ratio = _compute_fuzz_ratio(working_row['Metabolite'], row['Metabolite'])
        if ratio:
            ratios[row['Metabolite']] = ratio
    return ratios


def lipid_match(working_row, df):
    """
    """
    working_df = df[~df.index.isin([working_row.name])]
    wrow_is_lipid = not pandas.isna(working_row['lipid_nums'])
    matches = set()
    for index, row in working_df.iterrows():
        row_is_lipid = not pandas.isna(row['lipid_nums'])
        if wrow_is_lipid and row_is_lipid and working_row['lipid_name'] == row['lipid_name'] and working_row['lipid_nums'] == row['lipid_nums']:
            matches.add(row['Metabolite'])
    return matches


def isomer_match(working_row, df):
    """
    """
    working_df = df[~df.index.isin([working_row.name])]
    isomer_in_wrow = 'isomer' in working_row['lowered_name']
    wrow_starts_with_iso = working_row['lowered_name'].startswith('iso')
    wrow_starts_with_beta = working_row['lowered_name'].startswith('beta')
    matches = set()
    for index, row in working_df.iterrows():
        # MG(20:5)isomer and MG(20:5)
        isomer_in_row = 'isomer' in row['lowered_name']
        # Nicotinic acid and Isonicotinic acid
        row_starts_with_iso = row['lowered_name'].startswith('iso')
        # MALTOSE and BETA-MALTOSE
        row_starts_with_beta = row['lowered_name'].startswith('beta')
        lipid_matched = row['Metabolite'] in working_row['lipid_matches']
        fuzzy_matched = row['Metabolite'] in working_row['fuzzy_matches']
        in_matched = working_row['lowered_name'] in row['lowered_name'] or row['lowered_name'] in working_row['lowered_name']
        isomer_name_matched = row['isomer_name'] == working_row['isomer_name']
        both_have_isomer_nums = not (pandas.isna(row['isomer_nums']) and pandas.isna(working_row['isomer_nums']))
        isomer_nums_matched = both_have_isomer_nums and row['isomer_nums'] == working_row['isomer_nums']
        iso_counts_match = (not row['has_iso_count'] and not working_row['has_iso_count']) or (working_row['iso_counts'] == row['iso_counts'])
        matched = lipid_matched or fuzzy_matched or in_matched
        if iso_counts_match and \
           (((isomer_in_row or isomer_in_wrow) and matched) or \
            (isomer_name_matched and not isomer_nums_matched) or \
            ((wrow_starts_with_iso or row_starts_with_iso) and in_matched) or \
            ((wrow_starts_with_beta or row_starts_with_beta) and in_matched)\
           ):
            matches.add(row['Metabolite'])
    return matches


def _isotopologue_match(working_row, df, max_len_diff):
    """
    """
    working_df = df[~df.index.isin([working_row.name])]
    matches = set()
    # has_iso_count = any(count is not None for count in working_row['iso_counts'])
    # if has_iso_count:
    if working_row['has_iso_count']:
        working_len = len(working_row['iso_name'])
        for index, row in working_df.iterrows():
            iso_counts_match = working_row['iso_counts'] == row['iso_counts']
            name_in_name = working_row['iso_name'] in row['iso_name']
            names_fuzzy_match = row['fuzzy_matches'].get(working_row['Metabolite'])
            row_len = len(row['iso_name'])
            len_diff = abs(working_len - row_len)
            names_are_similar = (len_diff < max_len_diff and len_diff < min(row_len, working_len) and name_in_name) or names_fuzzy_match
            if iso_counts_match and names_are_similar:
                matches.add(row['Metabolite'])
    
    return matches if matches else ''

def _isotopologue_match_propagate(working_row, df):
    """
    """
    working_df = df[~df.index.isin([working_row.name])]
    matches = set() if not working_row['isotopologue_duplicates'] else working_row['isotopologue_duplicates']
    for index, row in working_df.iterrows():
        if working_row['Metabolite'] in row['isotopologue_duplicates']:
            matches.add(row['Metabolite'])
    
    return matches if matches else ''


# TODO add similar code for chemical formula, but it needs to parse and count the molecules in the formula, not just compare directly.
def _metadata_match(working_row, df, column_names_to_match, name_column='Metabolite'):
    """
    """
    working_df = df[~df.index.isin([working_row.name])]
    matches = {}
    for column_name in column_names_to_match:
        column = working_df.loc[:, column_name]
        match_mask = column.isin([working_row[column_name]])
        not_na_mask = ~column.isna()
        matching_indexes = column[(match_mask & not_na_mask)].index.values
        if matching_indexes.size > 0:
            matches[column_name] = working_df.loc[matching_indexes, name_column].values
    
    # TODO not sure if multiple mz/rt columns should be ORed or ANDed. Going with OR, for now, needs tested.
    return {match for column_name, matches in matches.items() for match in matches}



def _add_group_num(df, group_column_name, collection_column_name, metabolite_column_name='Metabolite'):
    """
    """
    name_to_group = {}
    group_num = 1
    df.loc[:, group_column_name] = pandas.Series([pandas.NA]*df.shape[0], index=df.index)
    for index, row in df.iterrows():
        name = row[metabolite_column_name]
        if dups := row[collection_column_name]:
            group_found = False
            for dup in dups:
                if dup in name_to_group:
                    df.loc[index, group_column_name] = name_to_group[dup]
                    name_to_group[name] = name_to_group[dup]
                    group_found = True
            if not group_found:
                df.loc[index, group_column_name] = group_num
                name_to_group[name] = group_num
                for dup in dups:
                    name_to_group[dup] = group_num
                    df.loc[df[df.loc[:, metabolite_column_name].isin([dup])].index, group_column_name] = group_num
                group_num += 1
    return df


def add_info_match_columns(df, met_dup_group_df):
    """
    """
    is_a_subset = pandas.Series(_determine_row_subsets(met_dup_group_df.drop(['Metabolite'], axis=1)), index=met_dup_group_df.index)
    no_subset_df = met_dup_group_df[~is_a_subset]
    for_info_match = no_subset_df.drop(['Metabolite'], axis=1) if no_subset_df.shape[1] > 1 else no_subset_df
    if no_subset_df.shape[1] > 1:
        info_match_grp = pandas.Series([set([num]) for num in pandas.factorize(for_info_match.apply(tuple, axis=1))[0]], 
                                       index = no_subset_df.index, 
                                       name='info_match_group')
    else:
        info_match_grp = pandas.Series([set([0]) for i in range(no_subset_df.shape[0])], index = no_subset_df.index, name='info_match_group')
    df = df.join(info_match_grp)
    no_info_group_mask = df.loc[:, 'info_match_group'].isna()
    no_info_grp_df = met_dup_group_df[no_info_group_mask]
    df.loc[no_info_group_mask, 'info_match_group'] = pandas.Series([set() for num in range(no_info_group_mask[no_info_group_mask].shape[0])],
                                                                                index = df.loc[no_info_group_mask, :].index)
    no_subset_row_sets = [set([(col, value) for col, value in row.items() if not pandas.isna(value)]) for index, row in for_info_match.iterrows()]
    for index, row in no_info_grp_df.drop(['Metabolite'], axis=1).iterrows():
        row_as_set = set([(col, value) for col, value in row.items() if not pandas.isna(value)])
        for j, row_set in enumerate(no_subset_row_sets):
            if row_as_set.issubset(row_set):
                df.loc[index, 'info_match_group'].update(df.loc[no_subset_df.index[j], 'info_match_group'])
    return df

def add_name_match_columns(df):
    """
    """
    df.loc[:, 'str_len'] = df.loc[:, 'Metabolite'].str.len()
    df.loc[:, 'lowered_name'] = df.loc[:, 'Metabolite'].str.lower().astype('string[pyarrow]')
    
    df.loc[:, 'iso_counts'] = df.loc[:, 'Metabolite'].apply(get_isotope_counts)
    df.loc[:, 'has_iso_count'] = df.apply(lambda row: any(count is not None for count in row['iso_counts']), axis=1, result_type='reduce').astype(bool)
    df.loc[:, 'iso_name'] = df.loc[:, 'lowered_name'].apply(remove_isotope_substring)
    
    df.loc[:, 'lipid_nums'] = df.loc[:, 'lowered_name'].str.extractall(r'([0-9a-z]+:[0-9a-z]+)').groupby(level=0)[0].apply('|'.join)
    df.loc[:, 'lipid_name'] = df.loc[:, 'lowered_name'].str.replace(r'(\(([^)]+:.+)\))|([0-9a-z]+:[0-9a-z]+[/_-])|([0-9a-z]+:[0-9a-z]+)', '', regex=True)
    df.loc[:, 'lipid_name'] = df.loc[:, 'lipid_name'].str.replace(r'(@\d+\.\d+$)|([ _-]\d$)|(\(\d\)$)|(\(.*isomer.*\)$)|(isomer\d*$)', '', regex=True)
    df.loc[:, 'lipid_name'] = df.loc[:, 'lipid_name'].str.strip()
    df.loc[:, 'lipid_matches'] = df.apply(lipid_match, args=(df,), axis=1, result_type='reduce').astype(object)
    
    # TODO remove the magic number 5 from this line, give it a decent name.
    df.loc[:, 'fuzzy_matches'] = df.apply(fuzzy_match, args=(df, 5), axis=1, result_type='reduce').astype(object)
    
    df.loc[:, 'CAS#'] = df.loc[:, 'lowered_name'].str.extract(r'\(cas# *(\d+-\d+-\d+) *\)')
    df.loc[:, 'short_name'] = df.loc[:, 'lipid_name'].apply(remove_isotope_substring).astype('string[pyarrow]')
    df.loc[:, 'short_name'] = df.loc[:, 'short_name'].str.replace(r'\(cas# *(\d+-\d+-\d+) *\)', '', regex=True)
    df.loc[:, 'short_name'] = df.loc[:, 'short_name'].str.replace(r'\(((possibly)|(\?+))\)', '', regex=True)
    df.loc[:, 'short_name'] = df.loc[:, 'short_name'].str.strip()
    df.loc[:, 'isomer_nums'] = df.loc[:, 'short_name'].str.extractall(r"(?:[^a-z0-9]|^)((?:\d{1,2}['`^]?(?:[,_]\d{1,2}['`^]?)*) *)-").groupby(level=0)[0].apply('|'.join)
    df.loc[:, 'isomer_nums'] = df.loc[:, 'isomer_nums'].str.replace(r"'|`|\^", '', regex=True)
    df.loc[:, 'isomer_nums'] = df.loc[:, 'isomer_nums'].str.replace(r"_", ',', regex=True)
    # df.loc[:, 'bond_nums'] = df.loc[:, 'short_name'].str.extractall(r"(\d{1,2}[ez])").groupby(level=0)[0].apply('|'.join)
    df.loc[:, 'isomer_name'] = df.loc[:, 'short_name'].str.replace(r"(?:[^a-z0-9]|^)((\d{1,2}['`^]?(,\d{1,2}['`^]?)*) *-)", '', regex=True)
    df.loc[:, 'isomer_name'] = df.loc[:, 'isomer_name'].str.replace(r"(\d{1,2}[ez])", '', regex=True)
    df.loc[:, 'isomer_matches'] = df.apply(isomer_match, args=(df,), axis=1, result_type='reduce').astype(object)
    
    # TODO remove the magic number 4 from this line, give it a decent name.
    df.loc[:, 'isotopologue_duplicates'] = df.apply(_isotopologue_match, args=(df, 4), axis=1, result_type='reduce').astype(object)
    df.loc[:, 'isotopologue_duplicates'] = df.apply(_isotopologue_match_propagate, args=(df,), axis=1, result_type='reduce').astype(object)
    df = _add_group_num(df, 'iso_group_num', 'isotopologue_duplicates')
    df = _add_group_num(df, 'fuzzy_group_num', 'fuzzy_matches')
    df = _add_group_num(df, 'isomer_group_num', 'isomer_matches')
    df = _add_group_num(df, 'lipid_group_num', 'lipid_matches')
    return df

def add_mz_rt_match_columns(df, column_finders):
    """
    """
    column_name_map = {column:column.lower() for column in df.columns}
    mz_column_names = column_finders['moverz_quant'].name_dict_match(column_name_map)
    rt_column_names = column_finders['retention_time'].name_dict_match(column_name_map)
    df.loc[:, 'mz_matches'] = df.apply(_metadata_match, args=(df, mz_column_names), axis=1, result_type='reduce').astype(object)
    df.loc[:, 'rt_matches'] = df.apply(_metadata_match, args=(df, rt_column_names), axis=1, result_type='reduce').astype(object)
    df.loc[:, 'mz_NA'] = df.apply(lambda row: row[mz_column_names].isna().all(), axis=1, result_type='reduce').astype(bool)
    df.loc[:, 'rt_NA'] = df.apply(lambda row: row[rt_column_names].isna().all(), axis=1, result_type='reduce').astype(bool)
    df.loc[:, 'mz_rt_NA'] = (df.loc[:, 'mz_NA']) & (df.loc[:, 'rt_NA'])
    if ~df.loc[:, 'mz_NA'].all() and df.loc[:, 'rt_NA'].all():
        df.loc[:, 'mz_rt_matches'] = df.loc[:, 'mz_matches']
    elif df.loc[:, 'mz_NA'].all() and ~df.loc[:, 'rt_NA'].all():
        df.loc[:, 'mz_rt_matches'] = df.loc[:, 'rt_matches']
    else:
        df.loc[:, 'mz_rt_matches'] = df.apply(lambda row: row['mz_matches'].intersection(row['rt_matches']), axis=1, result_type='reduce').astype(object)
    df = _add_group_num(df, 'mz_rt_group_num', 'mz_rt_matches')

    return df

def add_formula_match_columns(df, column_finders):
    """
    """
    column_name_map = {column:column.lower() for column in df.columns}
    formula_column_names = column_finders['formula'].name_dict_match(column_name_map)
    df.loc[:, 'formula_matches'] = df.apply(_metadata_match, args=(df, formula_column_names), axis=1, result_type='reduce').astype(object)
    df.loc[:, 'formula_NA'] = df.apply(lambda row: row[formula_column_names].isna().all(), axis=1, result_type='reduce').astype(bool)
    df = _add_group_num(df, 'formula_group_num', 'formula_matches')

    return df


def _group_num_to_final_group(df, final_group_num, group_column_name, 
                              info_group_column_name = 'info_match_group', 
                              final_group_column_name = 'final_group_num'):
    """
    """
    group_nums = {num for num in df.loc[:, group_column_name].values if not pandas.isna(num)}
    group_nums = sorted(list(group_nums))
    for group_num in group_nums:
        group_mask = df.loc[:, group_column_name].isin([group_num])
        final_group_mask = df.loc[:, final_group_column_name].isna()
        group_df = df[(group_mask & final_group_mask)]
        info_group_nums = {num for info_set in group_df.loc[:, info_group_column_name].values for num in info_set}
        info_group_nums = sorted(list(info_group_nums))
        for info_group_num in info_group_nums:
            info_group_mask = group_df.loc[:, info_group_column_name].apply(lambda x: info_group_num in x)
            sub_final_group_mask = df.loc[group_df.index, final_group_column_name].isna()
            final_group_df = group_df[(info_group_mask & sub_final_group_mask)]
            if final_group_df.shape[0]:
                df.loc[final_group_df.index, final_group_column_name] = final_group_num
                final_group_num += 1
    return df

def add_final_group_num(df, group_start_num=1):
    """
    """
    final_group_num = group_start_num
    df.loc[:, 'final_group_num'] = pandas.Series([pandas.NA]*df.shape[0])
    
    df = _group_num_to_final_group(df, group_start_num, 'iso_group_num')
    final_group_num = df.loc[:, 'final_group_num'].max()+1
    final_group_num = final_group_num if not pandas.isna(final_group_num) else group_start_num
    df = _group_num_to_final_group(df, final_group_num, 'isomer_group_num')
    final_group_num = df.loc[:, 'final_group_num'].max()+1
    final_group_num = final_group_num if not pandas.isna(final_group_num) else group_start_num
    df = _group_num_to_final_group(df, final_group_num, 'fuzzy_group_num')
    final_group_num = df.loc[:, 'final_group_num'].max()+1
    final_group_num = final_group_num if not pandas.isna(final_group_num) else group_start_num
    df = _group_num_to_final_group(df, final_group_num, 'mz_rt_group_num')
    final_group_num = df.loc[:, 'final_group_num'].max()+1
    final_group_num = final_group_num if not pandas.isna(final_group_num) else group_start_num
    df = _group_num_to_final_group(df, final_group_num, 'formula_group_num')
    
    # Set the rest to their own groups.
    final_group_mask = df.loc[:, 'final_group_num'].isna()
    final_group_num = df.loc[:, 'final_group_num'].max()+1
    final_group_num = final_group_num if not pandas.isna(final_group_num) else group_start_num
    for index, row in df[final_group_mask].iterrows():
        df.loc[index, 'final_group_num'] = final_group_num
        final_group_num += 1
    
    return df

def add_duplicate_rank(df, sort_by):
    """
    """
    df.loc[:, 'duplicate_ranks'] = pandas.Series(['']*df.shape[0], index=df.index, dtype=object)
    final_group_nums = {num for num in df.loc[:, 'final_group_num'].values}
    final_group_nums = sorted(list(final_group_nums))
    for final_group_num in final_group_nums:
        group_df = df[df.loc[:, 'final_group_num'].isin([final_group_num])]
        group_df = group_df.sort_values(by=sort_by)
        group_df.loc[:, 'duplicate_ranks'] = [str(num) for num in range(1, group_df.shape[0]+1)]
        df.update(group_df.loc[:, 'duplicate_ranks'])
    return df

def add_duplicate_attributes(df):
    """
    """
    df.loc[:, 'duplicate_attributes'] = pandas.Series([{'type':'Unknown'} for i in range(df.shape[0])], index=df.index, dtype=object)
    final_groups = df.groupby('final_group_num', dropna=False)
    for name, final_group in final_groups:
        if len(final_group) == 1:
            continue
        has_iso_group = ~final_group['iso_group_num'].isna().all()
        has_iso_count = final_group['has_iso_count'].any()
        has_fuzzy_group = ~final_group['fuzzy_group_num'].isna().all()
        has_mz_rt_group = ~final_group['mz_rt_group_num'].isna().all()
        has_isomer_group = ~final_group['isomer_group_num'].isna().all()
        mz_rt_NA = final_group['mz_rt_NA'].all()
        
        # If the lipid groups aren't all NA or the same group, but lipid nums exist, then the group cannot be equivalent.
        # same_lipid_nums = final_group['lipid_nums'].unique().shape[0] == 1
        same_lipid_group = final_group['lipid_group_num'].unique().shape[0] == 1
        lipid_nums_NA = final_group['lipid_nums'].isna().all()
        lipid_group_NA = final_group['lipid_group_num'].isna().all()
        equivalent_lipid = (lipid_nums_NA and lipid_group_NA) or (not lipid_group_NA and same_lipid_group)
        
        # If the isomer groups aren't all NA or the same group, but isomer nums exist, then the group cannot be equivalent.
        # same_isomer_nums = final_group['isomer_nums'].unique().shape[0] == 1
        same_isomer_group = final_group['isomer_group_num'].unique().shape[0] == 1
        isomer_nums_NA = final_group['isomer_nums'].isna().all()
        isomer_group_NA = final_group['isomer_group_num'].isna().all()
        equivalent_isomer = (isomer_nums_NA and isomer_group_NA) or (not isomer_group_NA and same_isomer_group)
                
        if equivalent_lipid and \
           equivalent_isomer and \
           (has_iso_group or (has_fuzzy_group and not has_iso_count)) and (mz_rt_NA or has_mz_rt_group):
            df.loc[final_group.index, 'duplicate_attributes'] = [{'type': 'Equivalent'} for i in range(final_group.shape[0])]
        elif equivalent_lipid and \
           (has_isomer_group or (has_mz_rt_group and not (has_iso_group or has_fuzzy_group))):
            df.loc[final_group.index, 'duplicate_attributes'] = [{'type': 'Ambiguous'} for i in range(final_group.shape[0])]
    return df


def _get_relevant_group_names(group, data_df, met_to_root):
    """
    """
    names = list(data_df.loc[group.index, :].iloc[:, 0])
    # Filter out names that are family members like a-glycerol-3P_2-13C0 and a-glycerol-3P_2-ALL   AN000186
    roots_to_names = {}
    for name in names:
        root = met_to_root[name]
        if root in roots_to_names:
            roots_to_names[root].add(name)
        else:
            roots_to_names[root] = set([name])
    names = [list(names)[0] for root, names in roots_to_names.items() if len(names) == 1]
    return names


"""
Make sure to note in the documentation that duplicates are found in DATA, and information in METABOLITES
is based on that. So a metabolite in METABOLIES marked as having the same information as another one 
is only because those metabolites have duplicate sample data. There could be additional metabolites 
that also have the same column information that are not marked because the sample data is not duplicated.
"""

# met_df = pandas.DataFrame({'Metabolite': ['alanine-13C3-15N', 'L-ALANINE_13C3', 'alanine-13c3-15n0', 'r-alanine', 'r_alanine', 'alanine', 'not-alanine', 'NOT-alanine', 'also_not_alanine', 'first_ambiguous', 'second_ambiguous'], 
#                            'col1': ['1', '1', pandas.NA, '1', pandas.NA, pandas.NA, '2', '2', '1', pandas.NA, pandas.NA],
#                            'col2': ['2', '2', pandas.NA, pandas.NA, pandas.NA, pandas.NA, '3', '3', '4', pandas.NA, pandas.NA],
#                            'mz': [pandas.NA, pandas.NA, pandas.NA, pandas.NA, pandas.NA, pandas.NA, pandas.NA, pandas.NA, '1', '10', '10'],
#                            'rt': [pandas.NA, pandas.NA, pandas.NA, pandas.NA, pandas.NA, pandas.NA, pandas.NA, pandas.NA, '4', '10', '10']})

# data_df = pandas.DataFrame({'Metabolite': ['alanine-13C3-15N', 'L-ALANINE_13C3', 'alanine-13c3-15n0', 'r-alanine', 'r_alanine', 'alanine', 'not-alanine', 'NOT-alanine', 'not_in_met', 'also_not_alanine', 'first_ambiguous', 'second_ambiguous'], 
#                            'col1': ['1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1'],
#                            'col2': ['2', '2', '2', '2', '2', '2', '2', '2', '2', '2', '2', '2']})

# met_df = pandas.DataFrame({'Metabolite': ["2'-asdf-4'-qwer", "2'-asdf-5'-qwer", "2'-asdf-6'-qwer", "2'-asdf-6'-qwer 13C1"], 
#                             'mz': ['1', '1', pandas.NA, '3']})

# data_df = pandas.DataFrame({'Metabolite': ["2'-asdf-4'-qwer", "2'-asdf-5'-qwer", "2'-asdf-6'-qwer", "2'-asdf-6'-qwer 13C1"], 
#                             'col1': ['1', '1', '1', '1'],
#                             'col2': ['2', '2', '2', '2']})


def augment_duplicates(data_df, met_df):
    """
    """
    data_duplicates = get_duplicate_rows(data_df, 'DATA')
    met_to_root, _, _ = find_metabolite_families(list(data_df.loc[:, "Metabolite"]))
    data_duplicate_mets = data_df.loc[data_duplicates.index, 'Metabolite']
    met_rows_of_data_dups = met_df[met_df.loc[:, 'Metabolite'].isin(data_duplicate_mets)]
    groups = data_duplicates.groupby(data_duplicates.columns.tolist(), dropna=False)
    names_in_data_not_met = []
    # Store data frames for each group to update the main dataframes with.
    data_group_dfs = []
    met_group_dfs = []
    for i, (name, group) in enumerate(groups):
        names = _get_relevant_group_names(group, data_df, met_to_root)
        if not names:
            continue
        
        data_dup_group_df = data_df.loc[group.index, :].loc[:, ['Metabolite']]
        data_dup_group_df = data_dup_group_df[data_dup_group_df.loc[:, 'Metabolite'].isin(names)]
        
        met_dup_group_df = met_rows_of_data_dups[met_rows_of_data_dups.loc[:, 'Metabolite'].isin(names)]
        intermediate_df = met_dup_group_df.copy()
        
        # If the names aren't in the METABOLITES table remove them from the DATA table.
        if len(met_dup_group_df.loc[:, 'Metabolite'].unique()) < len(set(names)):
            # Note, if somehow all the names aren't in METABOLITES, this will add them all to the list.
            names_in_data_not_met += [name for name in names if name not in met_dup_group_df.loc[:, "Metabolite"].values]
        
        
        intermediate_df = add_info_match_columns(intermediate_df, met_dup_group_df)
        intermediate_df = add_name_match_columns(intermediate_df)
        intermediate_df = add_mz_rt_match_columns(intermediate_df, column_finders)
        intermediate_df = add_formula_match_columns(intermediate_df, column_finders)
        intermediate_df.loc[:, 'num_NA_values'] = intermediate_df.isna().sum(axis=1)
        intermediate_df = add_final_group_num(intermediate_df, group_start_num=1)    
        intermediate_df = add_duplicate_rank(intermediate_df, ['num_NA_values', 'str_len', 'Metabolite'])
        intermediate_df = add_duplicate_attributes(intermediate_df)
        # intermediate_df.to_csv('C:/Users/Sparda/Desktop/Moseley Lab/Code/mwtab/scratch/temp1.csv', index=False)
            
        # Create the DATA intermediate_df.
        # If there are Metabolites in the METABOLITES section with the same name, we have to choose 1 to link back to DATA.
        # There has to be a 1-1 relationship. We can't add rows to DATA.
        sorted_intermediate_df = intermediate_df.sort_values(['Metabolite', 'duplicate_ranks'])
        dup_mets_removed = sorted_intermediate_df.loc[:, 'Metabolite'].drop_duplicates()
        intermediate_df_no_dups = intermediate_df.loc[dup_mets_removed.index, :]
        data_intermediate_df = data_dup_group_df.merge(intermediate_df_no_dups, on='Metabolite', how='left')
        data_intermediate_df.index = data_dup_group_df.index
        names_in_data_not_met_mask = data_dup_group_df.loc[:, 'Metabolite'].isin(names_in_data_not_met)
        if names_in_data_not_met_mask.any():
            data_not_in_met_df = data_dup_group_df.copy()
            data_not_in_met_df = data_not_in_met_df.loc[names_in_data_not_met_mask, :]
            data_not_in_met_df = add_name_match_columns(data_not_in_met_df)
            data_not_in_met_df = add_mz_rt_match_columns(data_not_in_met_df, column_finders)
            data_not_in_met_df = add_formula_match_columns(data_not_in_met_df, column_finders)
            final_group_num = intermediate_df.loc[:, 'final_group_num'].max()+1
            final_group_num = final_group_num if not pandas.isna(final_group_num) else 0
            data_not_in_met_df = add_final_group_num(data_not_in_met_df, final_group_num)
            data_not_in_met_df = add_duplicate_rank(data_not_in_met_df, ['str_len', 'Metabolite'])
            data_not_in_met_df = add_duplicate_attributes(data_not_in_met_df)
            data_not_in_met_df.loc[:, 'num_NA_values'] = pandas.NA
            data_intermediate_df = data_intermediate_df.loc[:, data_not_in_met_df.columns]
            data_intermediate_df.update(data_not_in_met_df)
        
        data_df_to_add = pandas.DataFrame({"duplicate_groups": [i+1]*data_intermediate_df.shape[0],
                                           "duplicate_sub_groups": data_intermediate_df.loc[:, 'final_group_num'].astype(int).values,
                                           "duplicate_ranks": data_intermediate_df.loc[:, 'duplicate_ranks'].values,
                                           "duplicate_attributes": data_intermediate_df.loc[:, 'duplicate_attributes'].astype(str).values},
                                          index = data_intermediate_df.index, dtype=str)
        data_group_dfs.append(data_df_to_add)
        
        met_df_to_add = pandas.DataFrame({"duplicate_groups": [i+1]*intermediate_df.shape[0],
                                          "duplicate_sub_groups": intermediate_df.loc[:, 'final_group_num'].values,
                                          "duplicate_ranks": intermediate_df.loc[:, 'duplicate_ranks'].values,
                                          "duplicate_attributes": intermediate_df.loc[:, 'duplicate_attributes'].astype(str).values},
                                         index = intermediate_df.index, dtype=str)
        met_group_dfs.append(met_df_to_add)
            
    
    if data_group_dfs:
        # TODO if Sub_Groups are identical to Ranks or duplicate_groups, then drop the column.
        data_df.loc[:, 'duplicate_groups'] = ''
        data_df.loc[:, 'duplicate_sub_groups'] = ''
        data_df.loc[:, 'duplicate_ranks'] = ''
        data_df.loc[:, 'duplicate_attributes'] = ''
        for df in data_group_dfs:
            data_df.update(df)
        
    if met_group_dfs:
        met_df.loc[:, 'duplicate_groups'] = ''
        met_df.loc[:, 'duplicate_sub_groups'] = ''
        met_df.loc[:, 'duplicate_ranks'] = ''
        met_df.loc[:, 'duplicate_attributes'] = ''
        for df in met_group_dfs:
            met_df.update(df)
    
    return data_df, met_df
    
# data_df.to_csv('C:/Users/Sparda/Desktop/Moseley Lab/Code/mwtab/scratch/temp1.csv', index=False)
# met_df.to_csv('C:/Users/Sparda/Desktop/Moseley Lab/Code/mwtab/scratch/temp2.csv', index=False)


logger = logging.getLogger('asdf')
logger.setLevel(logging.INFO)

file_handler = logging.FileHandler("C:/Users/Sparda/Desktop/Moseley Lab/Code/mwtab/scratch/augment.log", mode="w", encoding="utf-8")
file_handler.setLevel(logging.INFO)
logger.addHandler(file_handler)

logger.write = lambda msg: logger.info(msg) if msg != '\n' else None


# augment duplicates for relevant AN-IDs
with open('C:/Users/Sparda/Desktop/Moseley Lab/Code/mwtab/scratch/AN_IDs_with_duplicates.json','r') as jsonFile:
    AN_IDs_with_duplicates = json.load(jsonFile)

print(datetime.datetime.now())
with redirect_stdout(logger), redirect_stderr(logger):
    text_to_save = ''
    dfs_to_save = {}
    for an_id in AN_IDs_with_duplicates:
        # an_id = 'AN005077'
        try:
            data_df = pandas.read_csv(f'C:/Users/Sparda/Desktop/Moseley Lab/Code/mwtab/scratch/data_sections/{an_id}_data.csv', index_col=0, dtype=str)
            met_df = pandas.read_csv(f'C:/Users/Sparda/Desktop/Moseley Lab/Code/mwtab/scratch/metabolites_sections/{an_id}_metabolites.csv', index_col=0, dtype=str)
        except:
            continue
        
        # print(f'ANID: {an_id}\n')
        data_duplicates = get_duplicate_rows(data_df, 'DATA')
        data_df, met_df = augment_duplicates(data_df, met_df)
        text_to_save += an_id + '\n'
        groups = data_duplicates.groupby(data_duplicates.columns.tolist(), dropna=False)
        for name, group in groups:
            names = list(data_df.loc[group.index, :].iloc[:, 0])
            text_to_save += repr(names) + '\n'
            text_to_save += met_df[met_df.loc[:, 'Metabolite'].isin(names)].to_string() + '\n\n'
        dfs_to_save[an_id] =  met_df[met_df.loc[:, 'duplicate_attributes'] != '']

with open('C:/Users/Sparda/Desktop/Moseley Lab/Code/mwtab/scratch/augmented_duplicates.txt', 'wb') as outFile:
    outFile.write(text_to_save.encode("utf-8"))
print(datetime.datetime.now())

all_metabolite_names = []
for key, value in dfs_to_save.items():
    all_metabolite_names += list(value.loc[:, 'Metabolite'])
all_metabolite_names = list(set(all_metabolite_names))
with open('C:/Users/Sparda/Desktop/Moseley Lab/Code/mwtab/scratch/_all_duplicated_metabolites.json', 'w') as jsonFile:
    jsonFile.write(json.dumps(all_metabolite_names, indent=2))


text_to_save = ''
for key, value in dfs_to_save.items():
    if value.loc[:, 'duplicate_attributes'].str.contains('Equivalent').any():
        text_to_save += key + '\n'
        value = value.sort_values(['duplicate_groups', 'duplicate_sub_groups', 'duplicate_ranks'])
        value = value[value.loc[:, 'duplicate_attributes'].str.contains('Equivalent')]
        text_to_save += value.to_string() + '\n\n'

with open('C:/Users/Sparda/Desktop/Moseley Lab/Code/mwtab/scratch/_equivalent_duplicates.txt', 'wb') as outFile:
    outFile.write(text_to_save.encode("utf-8"))


text_to_save = ''
for key, value in dfs_to_save.items():
    if value.loc[:, 'duplicate_attributes'].str.contains('Equivalent').any():
        text_to_save += key + '\n'
        value = value.sort_values(['duplicate_groups', 'duplicate_sub_groups', 'duplicate_ranks'])
        value = value[value.loc[:, 'duplicate_attributes'].str.contains('Ambiguous')]
        text_to_save += value.to_string() + '\n\n'

with open('C:/Users/Sparda/Desktop/Moseley Lab/Code/mwtab/scratch/_ambiguous_duplicates.txt', 'wb') as outFile:
    outFile.write(text_to_save.encode("utf-8"))


text_to_save = ''
for key, value in dfs_to_save.items():
    if value.shape[1] == 5:
        text_to_save += key + '\n'
        value = value.sort_values(['duplicate_groups', 'duplicate_sub_groups', 'duplicate_ranks'])
        text_to_save += value.to_string() + '\n\n'

with open('C:/Users/Sparda/Desktop/Moseley Lab/Code/mwtab/scratch/_no_info_duplicates.txt', 'wb') as outFile:
    outFile.write(text_to_save.encode("utf-8"))

    

no_inf_dfs = {key:value.copy() for key, value in dfs_to_save.items() if value.shape[1] == 5}
no_inf_dfs2 = {key:value.copy() for key, value in dfs_to_save.items() if value.shape[1] == 5}

text_to_save = ''
text_to_save2 = ''
for key, value in no_inf_dfs.items():
    if not value.equals(no_inf_dfs2[key]):
        text_to_save += key + '\n'
        value2 = no_inf_dfs2[key].sort_values(['duplicate_groups', 'duplicate_sub_groups', 'duplicate_ranks'])
        text_to_save += value2.to_string() + '\n\n'
    else:
        text_to_save2 += key + '\n'
        value2 = no_inf_dfs2[key].sort_values(['duplicate_groups', 'duplicate_sub_groups', 'duplicate_ranks'])
        text_to_save2 += value2.to_string() + '\n\n'

with open('C:/Users/Sparda/Desktop/Moseley Lab/Code/mwtab/scratch/_no_info_duplicates.txt', 'wb') as outFile:
    outFile.write(text_to_save.encode("utf-8"))
with open('C:/Users/Sparda/Desktop/Moseley Lab/Code/mwtab/scratch/_no_info_duplicates2.txt', 'wb') as outFile:
    outFile.write(text_to_save2.encode("utf-8"))


# TODO Α-Linolenic Acid(C18:3N3)   γ-Linolenic Acid(C18:3N6)   Need to add code to find \d+:\d+ in names and make sure they are equal to fuzzy match.
# Most are surrounded by paranthesis, but some aren't.
# Similar examples  (±)12-HEPE [(±)-12-hydroxy-5Z,8Z,10E,14Z,17Z-eicosapentaenoic acid]   (±)15-HEPE [(±)-15-hydroxy-5Z,8Z,11Z,13E,17Z-eicosapentaenoic acid]
# PC 30:3; [M+H]+; GPCho(10:0/20:3(5Z,8Z,11Z))    PC 30:3; [M+H]+; GPCho(12:0/18:3(6Z,9Z,12Z))  comma seperated list of \d+(Z|E) 
# PE 28:1; [M+H]+; GPEtn(10:0/18:1(11E))
# PE 18:0-16:1
# ([cd]\d+:\d[ne][\do])  C18:3N6 d60:2EO  or (\d+:\d+[ep])     
# if matches \(.*:.*\) then pull the whole thing out, do the same, but surrounded by spaces or end of string as well
# ([0-9a-z]+:[0-9a-z]+)

# PysoPE 20:5(2n isomer1)   PysoPE 20:5  If name1 in name2 and 'isomer' in name2 then it's a match.
# N-Acetylglucosamine  and N-Acetylgalactosamine  fuzzy match when they shouldn't. Maybe look to see if the ends match, but the middle doesn't 
# and discount it? This is the only example like this so far.

# AN004165  isoleucine-13C and isoleucine-13C6 were matching as equivalent due to fuzzy match. They shouldn't. Add to tests.

# Test different ways to pull out lipid numbers.
temp = ['PI(P-22:4_12:1)_2194',
 'LPC 18:1/0:0',
 'LysoPC(20:2(11Z,14Z))',
 'LysoPC 20:2(2n isomer2)',
 '72.0579@6.181999:1']
df = pandas.DataFrame(temp, columns = ['Metabolite'])

df = pandas.DataFrame(all_metabolite_names, columns = ['Metabolite'])
df = df.loc[df.loc[:, 'Metabolite'].str.contains(':'), :]
df.loc[:, 'lowered_name'] = df.loc[:, 'Metabolite'].str.lower()
df.loc[:, 'lipid_nums1'] = df.loc[:, 'lowered_name'].str.extract('\((.+:.+)\)')
df.loc[:, 'lipid_nums1'] = df.loc[:, 'lipid_nums1'].str.extractall('([0-9a-z]+:[0-9a-z]+)').groupby(level=0)[0].apply(' '.join)
df.loc[:, 'lipid_nums2'] = df.loc[:, 'lowered_name'].str.extract('([0-9a-z]+:[0-9a-z/-_]+)')
df.loc[:, 'lipid_nums3'] = df.loc[:, 'lowered_name'].str.extractall('([0-9a-z]+:[0-9a-z]+)').groupby(level=0)[0].apply(' '.join)
df.loc[:, 'lipid_nums4'] = df.loc[:, 'lowered_name'].str.extractall('(?:^|[ (_/-])([0-9a-z]+:[0-9a-z]+)').groupby(level=0)[0].apply(' '.join)
different_mask = (df.loc[:, 'lipid_nums3'] == df.loc[:, 'lipid_nums4']) | ((df.loc[:, 'lipid_nums3'].isna()) & (df.loc[:, 'lipid_nums4'].isna()))
df.loc[~different_mask, :]


# Test fuzzy matching based on lipid nums matching.
def fuzzy_match2(working_row, df, max_len_diff):
    """
    """
    working_df = df[~df.index.isin([working_row.name])]
    working_len = len(working_row['lowered_name'])
    ratios = {}
    for index, row in working_df.iterrows():
        row_len = len(row['lowered_name'])
        len_diff = abs( working_len - row_len)
        if working_row['lowered_name'] in row['lowered_name'] and len_diff < max_len_diff and len_diff < min(row_len, working_len):
            ratio = 100
        else:
            ratio = _compute_fuzz_ratio(working_row['Metabolite'], row['Metabolite'])
        if ratio:
            ratios[row['Metabolite']] = ratio
    return ratios


lipid_dfs = {key:value.copy() for key, value in dfs_to_save.items() if value.loc[:, 'Metabolite'].str.contains(':').any()}
text_to_save = ''
for an_id, df in lipid_dfs.items():
    df = df.loc[df.loc[:, 'Metabolite'].str.contains(':'), :]
    df.loc[:, 'lowered_name'] = df.loc[:, 'Metabolite'].str.lower().astype('string[pyarrow]')
    df.loc[:, 'lipid_nums'] = df.loc[:, 'lowered_name'].str.extractall('([0-9a-z]+:[0-9a-z]+)').groupby(level=0)[0].apply('|'.join)
    text_to_save += an_id + '\n'
    for name, group in df.groupby('duplicate_groups'):
        group.loc[:, 'fuzzy_matches'] = group.apply(fuzzy_match, args=(group, 5), axis=1, result_type='reduce').astype(object)
        group.loc[:, 'fuzzy_matches2'] = group.apply(fuzzy_match2, args=(group, 5), axis=1, result_type='reduce').astype(object)
        fuzzy_mask = group.apply(lambda row: True if row['fuzzy_matches'].keys() == row['fuzzy_matches2'].keys() else False, axis=1)
        if not fuzzy_mask.all():
            group = group[~fuzzy_mask]
            group = group.sort_values(['duplicate_groups', 'duplicate_sub_groups', 'duplicate_ranks'])
            text_to_save += group.to_string() + '\n\n'
    
    
    # df = df.loc[:, ['Metabolite', 'lipid_nums', 'duplicate_groups', 'duplicate_sub_groups']]
    # text_to_save += an_id + '\n'
    # df = df.sort_values(['duplicate_groups', 'duplicate_sub_groups'])
    # text_to_save += df.to_string() + '\n\n'
with open('C:/Users/Sparda/Desktop/Moseley Lab/Code/mwtab/scratch/_lipid_nums_fuzzy_match.txt', 'wb') as outFile:
    outFile.write(text_to_save.encode("utf-8"))


# Test removing lipid nums from lipid names.
lipid_dfs = {key:value.copy() for key, value in dfs_to_save.items() if value.loc[:, 'Metabolite'].str.contains(':').any()}
text_to_save = ''
for an_id, df in lipid_dfs.items():
    df = df.loc[df.loc[:, 'Metabolite'].str.contains(':'), :]
    df.loc[:, 'lowered_name'] = df.loc[:, 'Metabolite'].str.lower().astype('string[pyarrow]')
    df.loc[:, 'lipid_nums'] = df.loc[:, 'lowered_name'].str.extractall('([0-9a-z]+:[0-9a-z]+)').groupby(level=0)[0].apply('|'.join)
    df.loc[:, 'lipid_name'] = df.loc[:, 'lowered_name'].str.replace('(\(([^)]+:.+)\))|([0-9a-z]+:[0-9a-z]+[/_-])|([0-9a-z]+:[0-9a-z]+)', '', regex=True)
    df.loc[:, 'lipid_name'] = df.loc[:, 'lipid_name'].str.replace('(@\d+\.\d+$)|([ _-]\d$)|(\(\d\)$)|(\(.*isomer.*\)$)|(isomer\d*$)', '', regex=True)
    df.loc[:, 'lipid_name'] = df.loc[:, 'lipid_name'].str.strip()
    text_to_save += an_id + '\n'
    for name, group in df.groupby('duplicate_groups'):
        if group.loc[:, 'lipid_nums'].duplicated().any():
            group = group.sort_values(['duplicate_groups', 'duplicate_sub_groups', 'duplicate_ranks'])
            group = group.loc[:, ['Metabolite', 'lipid_nums', 'lipid_name']]
            text_to_save += group.to_string() + '\n\n'

with open('C:/Users/Sparda/Desktop/Moseley Lab/Code/mwtab/scratch/_lipid_names_test.txt', 'wb') as outFile:
    outFile.write(text_to_save.encode("utf-8"))


# Test removing isomer nums from names.
def name_in_other_names(working_row, df, column_name):
    working_df = df[~df.index.isin([working_row.name])]
    return any(working_row[column_name] in value for value in working_df.loc[:, column_name].values)

def values_equal(working_row, df, column_name):
    working_df = df[~df.index.isin([working_row.name])]
    return any(not pandas.isna(value) and not pandas.isna(working_row[column_name]) and working_row[column_name] == value for value in working_df.loc[:, column_name].values)


dfs = {key:value.copy() for key, value in dfs_to_save.items()}
text_to_save = ''
for an_id, df in dfs.items():
    df.loc[:, 'lowered_name'] = df.loc[:, 'Metabolite'].str.lower().astype('string[pyarrow]')
    df.loc[:, 'lipid_name'] = df.loc[:, 'lowered_name'].str.replace(r'(\(([^)]+:.+)\))|([0-9a-z]+:[0-9a-z]+[/_-])|([0-9a-z]+:[0-9a-z]+)', '', regex=True)
    df.loc[:, 'lipid_name'] = df.loc[:, 'lipid_name'].str.replace(r'(@\d+\.\d+$)|([ _-]\d$)|(\(\d\)$)|(\(.*isomer.*\)$)|(isomer\d*$)', '', regex=True)
    df.loc[:, 'lipid_name'] = df.loc[:, 'lipid_name'].str.strip()
    
    df.loc[:, 'CAS#'] = df.loc[:, 'lowered_name'].str.extract(r'\(cas# *(\d+-\d+-\d+) *\)')
    df.loc[:, 'short_name'] = df.loc[:, 'lipid_name'].apply(remove_isotope_substring).astype('string[pyarrow]')
    df.loc[:, 'short_name'] = df.loc[:, 'short_name'].str.replace(r'\(cas# *(\d+-\d+-\d+) *\)', '', regex=True)
    df.loc[:, 'short_name'] = df.loc[:, 'short_name'].str.replace(r'\(((possibly)|(\?+))\)', '', regex=True)
    df.loc[:, 'short_name'] = df.loc[:, 'short_name'].str.strip()
    df.loc[:, 'isomer_nums'] = df.loc[:, 'short_name'].str.extractall(r"(?:[^a-z0-9]|^)((?:\d{1,2}['`^]?(?:[,_]\d{1,2}['`^]?)*) *)-").groupby(level=0)[0].apply('|'.join)
    # df.loc[:, 'isomer_nums1'] = df.loc[:, 'short_name'].str.extractall(r"(?:[^a-z0-9]|^)((?:\d{1,2}['`^]?(?:[,_]\d{1,2}['`^]?)*) *)-").groupby(level=0)[0].apply('|'.join)
    # isomer_equal_mask = (df.loc[:, 'isomer_nums'] == df.loc[:, 'isomer_nums1']) | (df.loc[:, 'isomer_nums'].isna() & df.loc[:, 'isomer_nums1'].isna())
    # if not isomer_equal_mask.all():
    #     print(an_id)
    #     print(df.loc[~isomer_equal_mask, ['Metabolite', 'isomer_nums', 'isomer_nums1']])
    #     print()
    # Not sure if this is needed yet, asking Hunter if those characters after numbers mean anything.
    df.loc[:, 'isomer_nums'] = df.loc[:, 'isomer_nums'].str.replace(r"'|`|\^", '', regex=True)
    df.loc[:, 'isomer_nums'] = df.loc[:, 'isomer_nums'].str.replace(r"_", ',', regex=True)
    df.loc[:, 'isomer_name'] = df.loc[:, 'short_name'].str.replace(r"(?:[^a-z0-9]|^)((\d{1,2}['`^]?(,\d{1,2}['`^]?)*) *-)", '', regex=True)
    # Need to check isomer_names in each other and isomer_nums equal. Also need to check isomer_name equal and isomer_nums different.
    text_to_save += an_id + '\n'
    for name, group in df.groupby('duplicate_groups'):
        isomer_name_in_isomer_name = group.apply(name_in_other_names, args=(group, 'isomer_name'), axis=1, result_type='reduce').any()
        isomer_nums_equal = group.apply(values_equal, args=(group, 'isomer_nums'), axis=1, result_type='reduce').any()
        isomer_name_equal = group.apply(values_equal, args=(group, 'isomer_name'), axis=1, result_type='reduce').any()
        names_in_other_names = group.apply(name_in_other_names, args=(group, 'Metabolite'), axis=1, result_type='reduce').any()
        if not group.loc[:, 'isomer_nums'].isna().all() and group.shape[0] > 1 and not names_in_other_names and isomer_name_equal and not isomer_nums_equal:
            group = group.sort_values(['duplicate_groups', 'duplicate_sub_groups', 'duplicate_ranks'])
            group = group.loc[:, ['Metabolite', 'isomer_nums', 'isomer_name']]
            text_to_save += group.to_string() + '\n\n'

with open('C:/Users/Sparda/Desktop/Moseley Lab/Code/mwtab/scratch/_isomer_names_test.txt', 'wb') as outFile:
    outFile.write(text_to_save.encode("utf-8"))

# TODO
# If isomer_name is equal and isomer_nums are equal, then mark as equivalent, unless lipid_nums, mz, rt, or formula contradict.
# If isomer_name is equal and isomer_nums are different, then mark as ambiguous, unless lipid_nums, mz, rt, or formula contradict.


# Testing remove bond_nums from name.
dfs = {key:value.copy() for key, value in dfs_to_save.items()}
text_to_save = ''
for an_id, df in dfs.items():
    df.loc[:, 'lowered_name'] = df.loc[:, 'Metabolite'].str.lower().astype('string[pyarrow]')
    df.loc[:, 'lipid_name'] = df.loc[:, 'lowered_name'].str.replace(r'(\(([^)]+:.+)\))|([0-9a-z]+:[0-9a-z]+[/_-])|([0-9a-z]+:[0-9a-z]+)', '', regex=True)
    df.loc[:, 'lipid_name'] = df.loc[:, 'lipid_name'].str.replace(r'(@\d+\.\d+$)|([ _-]\d$)|(\(\d\)$)|(\(.*isomer.*\)$)|(isomer\d*$)', '', regex=True)
    df.loc[:, 'lipid_name'] = df.loc[:, 'lipid_name'].str.strip()
    
    df.loc[:, 'short_name'] = df.loc[:, 'lipid_name'].apply(remove_isotope_substring).astype('string[pyarrow]')
    df.loc[:, 'short_name'] = df.loc[:, 'short_name'].str.replace(r'\(cas# *(\d+-\d+-\d+) *\)', '', regex=True)
    df.loc[:, 'short_name'] = df.loc[:, 'short_name'].str.replace(r'\(((possibly)|(\?+))\)', '', regex=True)
    df.loc[:, 'short_name'] = df.loc[:, 'short_name'].str.strip()
    df.loc[:, 'bond_nums'] = df.loc[:, 'short_name'].str.extractall(r"(\d{1,2}[ez])").groupby(level=0)[0].apply('|'.join)
    text_to_save += an_id + '\n'
    for name, group in df.groupby('duplicate_groups'):
        if not group.loc[:, 'bond_nums'].isna().all():
            group = group.sort_values(['duplicate_groups', 'duplicate_sub_groups', 'duplicate_ranks'])
            group = group.loc[:, ['Metabolite', 'bond_nums']]
            text_to_save += group.to_string() + '\n\n'

with open('C:/Users/Sparda/Desktop/Moseley Lab/Code/mwtab/scratch/_isomer_names_test.txt', 'wb') as outFile:
    outFile.write(text_to_save.encode("utf-8"))




