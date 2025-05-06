# -*- coding: utf-8 -*-
"""
Augment the DATA and METABOLITES tables of an mwtab file with extra columns.
"""
import re

import pandas

from . import utility_functions
from . import repair_metabolites_matching

_create_numeric_df = utility_functions._create_numeric_df
find_metabolite_families = utility_functions.find_metabolite_families
_determine_row_subsets = utility_functions._determine_row_subsets
compute_fuzz_ratios = utility_functions.compute_fuzz_ratios
_compute_fuzz_ratio = utility_functions._compute_fuzz_ratio

column_finders = repair_metabolites_matching.column_finders


def get_duplicate_rows(df, section_name, numeric=True):
    """
    """    
    # Convert every column except the metabolite column to numbers.
    if numeric:
        numeric_df = _create_numeric_df(df)
    else:
        numeric_df = df.iloc[:, 1:]
    
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
    compound = compound.lower()
    # re.sub(r'([^\d]|^|2h\d*|15n\d*)(?:(?:13c(?:_|\-)?(?:\d*))(15n|2h)|(?:13c(?:_|\-)?(?:\d*)))', _handle_replacement, 'cystine-13C615N2'.lower())
    # re.sub(r'([^\d]|^|2h\d*|15n\d*)(?:(?:13c(?:_|\-)?(?:\d*))(15n|2h)|(?:13c(?:_|\-)?(?:\d*)))', _handle_replacement, '13C-cystine'.lower())
    compound = re.sub(r'([^\d]|^|2h\d*|15n\d*)(?:(?:13c(?:_|\-)?(?:\d*))(15n|2h)|(?:13c(?:_|\-)?(?:\d*)))', _handle_replacement, compound)
    compound = re.sub(r'([^\d]|^|2h\d*|13c\d*)(?:(?:15n(?:_|\-)?(?:\d*))(13c|2h)|(?:15n(?:_|\-)?(?:\d*)))', _handle_replacement, compound)
    compound = re.sub(r'([^\d]|^|13c\d*|15n\d*)(?:(?:2h(?:_|\-)?(?:\d*))(15n|13c)|(?:2h(?:_|\-)?(?:\d*)))', _handle_replacement, compound)
    compound = compound.strip('_- ')
    return compound

def _handle_replacement(match_object):
    """
    """
    # print(match_object.groups())
    if match_object.group(2):
        return match_object.group(1) + match_object.group(2)
    return match_object.group(1)


def fuzzy_match(working_row, df):
    """
    """
    working_df = df[~df.index.isin([working_row.name])]
    ratios = {}
    for index, row in working_df.iterrows():
        ratio = _compute_fuzz_ratio(working_row['Metabolite'], row['Metabolite'])
        if ratio:
            ratios[row['Metabolite']] = ratio
    return ratios

# def fuzzy_match(iso_str, compare_list, short_name_to_long_name):
#     """
#     """
#     iso_str = iso_str.lower()
#     lowered_compare_list = [iso.lower() for iso in compare_list]
#     # lowered_to_original = {iso.lower():iso for iso in compare_list}
#     ratios = compute_fuzz_ratios([iso_str], [iso for iso in lowered_compare_list if iso != iso_str], inclusion_ratio=90)
#     if ratios:
#         fuzzy_series = ratios[iso_str]
#         fuzzy_dict = fuzzy_series.to_dict()
#         renamed_dict = {short_name_to_long_name[name]:ratio for name, ratio in fuzzy_dict.items()}
#     else:
#         return {}
#     return renamed_dict


def _isotopologue_match(working_row, df, max_len_diff):
    """
    """
    working_df = df[~df.index.isin([working_row.name])]
    matches = []
    has_iso_count = any(count is not None for count in working_row['iso_counts'])
    if has_iso_count:
        working_len = len(working_row['short_name'])
        for index, row in working_df.iterrows():
            iso_counts_match = working_row['iso_counts'] == row['iso_counts']
            name_in_name = working_row['short_name'] in row['short_name']
            names_fuzzy_match = row['fuzzy_matches'].get(working_row['Metabolite'])
            names_are_simliar_len = abs(working_len - len(row['short_name'])) < max_len_diff
            names_are_similar = ( names_are_simliar_len and name_in_name) or names_fuzzy_match
            if iso_counts_match and names_are_similar:
                matches.append(row['Metabolite'])
    
    return matches if matches else ''

def _isotopologue_match_propagate(working_row, df):
    """
    """
    working_df = df[~df.index.isin([working_row.name])]
    matches = [] if not working_row['isotopologue_duplicates'] else working_row['isotopologue_duplicates']
    for index, row in working_df.iterrows():
        if working_row['Metabolite'] in row['isotopologue_duplicates']:
            matches.append(row['Metabolite'])
    
    return matches if matches else ''



def _mz_rt_match(working_row, df, column_names_to_match, name_column='Metabolite'):
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
    # index_to_group = {}
    group_num = 0
    df.loc[:, group_column_name] = pandas.NA
    for index, row in df.iterrows():
        name = row[metabolite_column_name]
        if dups := row[collection_column_name]:
            group_found = False
            for dup in dups:
                if dup in name_to_group:
                    df.loc[index, group_column_name] = name_to_group[dup]
                    name_to_group[name] = name_to_group[dup]
                    # index_to_group[index] = name_to_group[dup]
                    group_found = True
            if not group_found:
                df.loc[index, group_column_name] = group_num
                name_to_group[name] = group_num
                # index_to_group[index] = group_num
                for dup in dups:
                    name_to_group[dup] = group_num
                    # dup_indexes = [dup_index for dup_index in df[df.loc[:, metabolite_column_name].isin([dup])].index]
                    # for index in dup_indexes:
                    #     index_to_group[index] = group_num
                    df.loc[df[df.loc[:, metabolite_column_name].isin([dup])].index, group_column_name] = group_num
                group_num += 1
    return df

def add_match_columns(df, column_finders):
    """
    """
    df.loc[:, 'str_len'] = df.loc[:, 'Metabolite'].str.len()
    df.loc[:, 'iso_counts'] = df.loc[:, 'Metabolite'].apply(get_isotope_counts)
    df.loc[:, 'short_name'] = df.loc[:, 'Metabolite'].apply(remove_isotope_substring)
    df.loc[:, 'fuzzy_matches'] = df.apply(fuzzy_match, args=(df,), axis=1)
    df.loc[:, 'isotopologue_duplicates'] = df.apply(_isotopologue_match, args=(df, 4), axis=1)
    df.loc[:, 'isotopologue_duplicates'] = df.apply(_isotopologue_match_propagate, args=(df,), axis=1)
    df = _add_group_num(df, 'iso_group_num', 'isotopologue_duplicates')
    df = _add_group_num(df, 'fuzzy_group_num', 'fuzzy_matches')
    
    column_name_map = {column:column.lower() for column in df.columns}
    mz_column_names = column_finders['moverz_quant'].name_dict_match(column_name_map)
    rt_column_names = column_finders['retention_time'].name_dict_match(column_name_map)
    df.loc[:, 'mz_matches'] = df.apply(_mz_rt_match, args=(df, mz_column_names), axis=1)
    df.loc[:, 'rt_matches'] = df.apply(_mz_rt_match, args=(df, rt_column_names), axis=1)
    df.loc[:, 'mz_rt_matches'] = df.apply(lambda row: row['mz_matches'].intersection(row['rt_matches']), axis=1)
    df = _add_group_num(df, 'mz_rt_group_num', 'mz_rt_matches')
    df.loc[:, 'mz_NA'] = df.apply(lambda row: row[mz_column_names].isna().all(), axis=1)
    df.loc[:, 'rt_NA'] = df.apply(lambda row: row[rt_column_names].isna().all(), axis=1)
    df.loc[:, 'mz_rt_NA'] = (df.loc[:, 'mz_NA']) & (df.loc[:, 'rt_NA'])
    
    
    
    # name_to_iso_group = {}
    # index_to_iso_group = {}
    # name_to_fuzzy_group = {}
    # index_to_fuzzy_group = {}
    # iso_group_num = 0
    # fuzzy_group_num = 0
    # df.loc[:, 'iso_group_num'] = GROUP_DEFAULT
    # df.loc[:, 'fuzzy_group_num'] = GROUP_DEFAULT
    # for index, row in df.iterrows():
    #     name = row['Metabolite']
    #     if iso_dups := row['isotopologue_duplicates']:
    #         iso_group_found = False
    #         for iso_dup in iso_dups:
    #             if iso_dup in name_to_iso_group:
    #                 df.loc[index, 'iso_group_num'] = name_to_iso_group[iso_dup]
    #                 name_to_iso_group[name] = name_to_iso_group[iso_dup]
    #                 index_to_iso_group[index] = name_to_iso_group[iso_dup]
    #                 iso_group_found = True
    #         if not iso_group_found:
    #             df.loc[index, 'iso_group_num'] = iso_group_num
    #             name_to_iso_group[name] = iso_group_num
    #             index_to_iso_group[index] = iso_group_num
    #             for iso_dup in iso_dups:
    #                 name_to_iso_group[iso_dup] = iso_group_num
    #                 iso_dup_indexes = [iso_index for iso_index in df[df.loc[:, 'Metabolite'].isin([iso_dup])].index]
    #                 for iso_index in iso_dup_indexes:
    #                     index_to_iso_group[iso_index] = iso_group_num
    #                 df.loc[df[df.loc[:, 'Metabolite'].isin([iso_dup])].index, 'iso_group_num'] = iso_group_num
    #             iso_group_num += 1
        
    #     if fuzzy_dups := row['fuzzy_matches']:
    #         fuzzy_group_found = False
    #         for fuzzy_dup in fuzzy_dups:
    #             if fuzzy_dup in name_to_fuzzy_group:
    #                 df.loc[index, 'fuzzy_group_num'] = name_to_fuzzy_group[fuzzy_dup]
    #                 name_to_fuzzy_group[name] = name_to_fuzzy_group[fuzzy_dup]
    #                 index_to_fuzzy_group[index] = name_to_fuzzy_group[fuzzy_dup]
    #                 fuzzy_group_found = True
    #         if not fuzzy_group_found:
    #             df.loc[index, 'fuzzy_group_num'] = fuzzy_group_num
    #             name_to_fuzzy_group[name] = fuzzy_group_num
    #             index_to_fuzzy_group[index] = fuzzy_group_num
    #             for fuzzy_dup in fuzzy_dups:
    #                 name_to_fuzzy_group[fuzzy_dup] = fuzzy_group_num
    #                 fuzzy_dup_indexes = [fuzzy_index for fuzzy_index in df[df.loc[:, 'Metabolite'].isin([fuzzy_dup])].index]
    #                 for fuzzy_index in fuzzy_dup_indexes:
    #                     index_to_fuzzy_group[fuzzy_index] = fuzzy_group_num
    #                 df.loc[df[df.loc[:, 'Metabolite'].isin([fuzzy_dup])].index, 'fuzzy_group_num'] = fuzzy_group_num
    #             fuzzy_group_num += 1
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

def add_final_group_num(df, group_start_num=0):
    """
    """
    final_group_num = group_start_num
    df.loc[:, 'final_group_num'] = pandas.NA
    
    df = _group_num_to_final_group(df, group_start_num, 'iso_group_num')
    final_group_num = df.loc[:, 'final_group_num'].max()+1
    final_group_num = final_group_num if not pandas.isna(final_group_num) else group_start_num
    df = _group_num_to_final_group(df, final_group_num, 'fuzzy_group_num')
    final_group_num = df.loc[:, 'final_group_num'].max()+1
    final_group_num = final_group_num if not pandas.isna(final_group_num) else group_start_num
    df = _group_num_to_final_group(df, final_group_num, 'mz_rt_group_num')
    
    # iso_group_nums = {num for num in df.loc[:, 'iso_group_num'].values if not pandas.isna(num)}
    # iso_group_nums = sorted(list(iso_group_nums))
    # for iso_group_num in iso_group_nums:
    #     iso_group_mask = df.loc[:, 'iso_group_num'].isin([iso_group_num])
    #     final_group_mask = df.loc[:, 'final_group_num'].isna()
    #     group_df = df[(iso_group_mask & final_group_mask)]
    #     info_group_nums = {num for info_set in group_df.loc[:, 'info_match_group'].values for num in info_set}
    #     info_group_nums = sorted(list(info_group_nums))
    #     for info_group_num in info_group_nums:
    #         info_group_mask = group_df.loc[:, 'info_match_group'].apply(lambda x: info_group_num in x)
    #         sub_final_group_mask = df.loc[group_df.index, 'final_group_num'].isna()
    #         final_group_df = group_df[(info_group_mask & sub_final_group_mask)]
    #         if final_group_df.shape[0]:
    #             df.loc[final_group_df.index, 'final_group_num'] = final_group_num
    #             final_group_num += 1
    
    # fuzzy_group_nums = {num for num in df.loc[:, 'fuzzy_group_num'].values if not pandas.isna(num)}
    # fuzzy_group_nums = sorted(list(fuzzy_group_nums))
    # for fuzzy_group_num in fuzzy_group_nums:
    #     fuzzy_group_mask = df.loc[:, 'fuzzy_group_num'].isin([fuzzy_group_num])
    #     final_group_mask = df.loc[:, 'final_group_num'].isna()
    #     group_df = df[(fuzzy_group_mask & final_group_mask)]
    #     info_group_nums = {num for info_set in group_df.loc[:, 'info_match_group'].values for num in info_set}
    #     info_group_nums = sorted(list(info_group_nums))
    #     for info_group_num in info_group_nums:
    #         info_group_mask = group_df.loc[:, 'info_match_group'].apply(lambda x: info_group_num in x)
    #         sub_final_group_mask = df.loc[group_df.index, 'final_group_num'].isna()
    #         final_group_df = group_df[(info_group_mask & sub_final_group_mask)]
    #         if final_group_df.shape[0]:
    #             df.loc[final_group_df.index, 'final_group_num'] = final_group_num
    #             final_group_num += 1
    
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
    df.loc[:, 'Duplicate Rank'] = ''
    final_group_nums = {num for num in df.loc[:, 'final_group_num'].values}
    final_group_nums = sorted(list(final_group_nums))
    for final_group_num in final_group_nums:
        group_df = df[df.loc[:, 'final_group_num'].isin([final_group_num])]
        group_df = group_df.sort_values(by=sort_by)
        group_df.loc[:, 'Duplicate Rank'] = [str(num) for num in range(group_df.shape[0])]
        df.update(group_df.loc[:, 'Duplicate Rank'])
    return df

def add_duplicate_attributes(df):
    """
    """
    df.loc[:, 'duplicate_attributes'] = [{'type':'Unknown'} for i in range(df.shape[0])]
    final_groups = df.groupby('final_group_num', dropna=False)
    for name, final_group in final_groups:
        has_iso_group = ~final_group['iso_group_num'].isna().all()
        has_fuzzy_group = ~final_group['fuzzy_group_num'].isna().all()
        has_mz_rt_group = ~final_group['mz_rt_group_num'].isna().all()
        mz_rt_NA = final_group['mz_rt_NA'].all()
        
        if (has_iso_group or has_fuzzy_group) and (mz_rt_NA or has_mz_rt_group):
            df.loc[final_group.index, 'duplicate_attributes'] = [{'type': 'Equivalent'} for i in range(final_group.shape[0])]
        if has_mz_rt_group and not (has_iso_group or has_fuzzy_group):
            df.loc[final_group.index, 'duplicate_attributes'] = [{'type': 'Ambiguous'} for i in range(final_group.shape[0])]
    return df


"""
I am thinking detecting duplicates should add 3 columns. 
1 like "duplicate_groups" that assigns a unique number to each group of duplicates.
1 like "Duplicates" that has the index and metabolite name of the others it is grouped with.
1 like "Duplicate_Info" that has the classification description or number for each duplicate.
The first 2 are redundant, I'm not sure both are needed.
Make sure to note in the documentation that duplicates are found in DATA, and information in METABOLITES
is based on that. So a metabolite in METABOLIES marked as having the same information as another one 
is only because those metabolites have duplicate sample data. There could be additional metabolites 
that also have the same column information that are not marked because the sample data is not duplicated.
"""
data_df = pandas.DataFrame()
met_df = pandas.DataFrame()

met_df = pandas.DataFrame({'Metabolite': ['alanine-13C3-15N', 'L-ALANINE_13C3', 'alanine-13c3-15n0', 'r-alanine', 'r_alanine', 'alanine', 'not-alanine', 'NOT-alanine', 'also_not_alanine', 'first_ambiguous', 'second_ambiguous'], 
                           'col1': ['1', '1', pandas.NA, '1', pandas.NA, pandas.NA, '2', '2', '1', pandas.NA, pandas.NA],
                           'col2': ['2', '2', pandas.NA, pandas.NA, pandas.NA, pandas.NA, '3', '3', '4', pandas.NA, pandas.NA],
                           'mz': [pandas.NA, pandas.NA, pandas.NA, pandas.NA, pandas.NA, pandas.NA, pandas.NA, pandas.NA, '1', '10', '10'],
                           'rt': [pandas.NA, pandas.NA, pandas.NA, pandas.NA, pandas.NA, pandas.NA, pandas.NA, pandas.NA, '4', '10', '10']})

data_df = pandas.DataFrame({'Metabolite': ['alanine-13C3-15N', 'L-ALANINE_13C3', 'alanine-13c3-15n0', 'r-alanine', 'r_alanine', 'alanine', 'not-alanine', 'NOT-alanine', 'not_in_met', 'also_not_alanine', 'first_ambiguous', 'second_ambiguous'], 
                           'col1': ['1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1'],
                           'col2': ['2', '2', '2', '2', '2', '2', '2', '2', '2', '2', '2', '2']})


# Remove duplicates from Data.
data_duplicates = get_duplicate_rows(data_df, 'DATA')
met_to_root, root_to_mets = find_metabolite_families(list(data_df.loc[:, "Metabolite"]))
data_duplicate_mets = data_df.loc[data_duplicates.index, 'Metabolite']
met_rows_of_data_dups = met_df[met_df.loc[:, 'Metabolite'].isin(data_duplicate_mets)]
groups = data_duplicates.groupby(data_duplicates.columns.tolist(), dropna=False)
names_in_data_not_met = []
# Store data frames for each group to update the main dataframes with.
data_group_dfs = []
met_group_dfs = []
for i, (name, group) in enumerate(groups):
    names = list(data_df.loc[group.index, :].iloc[:, 0])
    # Filter out names that are family members like a-glycerol-3P_2-13C0 and a-glycerol-3P_2-ALL   AN000186
    roots_to_names = {}
    for name2 in names:
        root = met_to_root[name2]
        if root in roots_to_names:
            roots_to_names[root].append(name2)
        else:
            roots_to_names[root] = [name2]
    names = [names[0] for root, names in roots_to_names.items() if len(names) == 1]
    original_names = names.copy()
    
    data_dup_group_df = data_df.loc[group.index, :].loc[:, ['Metabolite']]
    data_dup_group_df = data_dup_group_df[data_dup_group_df.loc[:, 'Metabolite'].isin(names)]
    
    met_dup_group_df = met_rows_of_data_dups[met_rows_of_data_dups.loc[:, 'Metabolite'].isin(names)]
    intermediate_df = met_dup_group_df.copy()
    
    # If the names aren't in the METABOLITES table remove them from the DATA table.
    if len(met_dup_group_df.loc[:, 'Metabolite'].unique()) < len(set(names)):
        # Note, if somehow all the names aren't in METABOLITES, this will add them all to the list.
        names_in_data_not_met += [name for name in names if name not in met_dup_group_df.loc[:, "Metabolite"].values]
    
    
    is_a_subset = pandas.Series(_determine_row_subsets(met_dup_group_df.iloc[:, 1:]))
    not_na_rows = ~met_dup_group_df.iloc[:, 1:].isna().all(axis=1)
    no_subset_df = met_dup_group_df[~is_a_subset]
    info_match_grp = pandas.Series([set([num]) for num in pandas.factorize(no_subset_df.iloc[:, 1:].apply(tuple, axis=1))[0]], 
                                   index = no_subset_df.index, 
                                   name='info_match_group')
    intermediate_df = intermediate_df.join(info_match_grp)
    no_info_group_mask = intermediate_df.loc[:, 'info_match_group'].isna()
    # no_info_grp_df = met_dup_group_df[(no_info_group_mask & not_na_rows)]
    no_info_grp_df = met_dup_group_df[no_info_group_mask]
    intermediate_df.loc[no_info_group_mask, 'info_match_group'] = pandas.Series([set() for num in range(no_info_group_mask[no_info_group_mask].shape[0])],
                                                                                index = intermediate_df.loc[no_info_group_mask, :].index)
    row_len = no_info_grp_df.shape[1] - 1
    no_subset_row_sets = [set([(col, value) for col, value in row.items() if not pandas.isna(value)]) for index, row in no_subset_df.iloc[:, 1:].iterrows()]
    for index, row in no_info_grp_df.iloc[:, 1:].iterrows():
        row_as_set = set([(col, value) for col, value in row.items() if not pandas.isna(value)])
        for j, row_set in enumerate(no_subset_row_sets):
            if row_as_set.issubset(row_set):
                intermediate_df.loc[index, 'info_match_group'].update(intermediate_df.loc[no_subset_df.index[j], 'info_match_group'])
    
    intermediate_df = add_match_columns(intermediate_df, column_finders)
    intermediate_df.loc[:, 'num_NA_values'] = intermediate_df.isna().sum(axis=1)
    intermediate_df = add_final_group_num(intermediate_df, group_start_num=0)    
    intermediate_df = add_duplicate_rank(intermediate_df, ['num_NA_values', 'str_len', 'Metabolite'])
    intermediate_df = add_duplicate_attributes(intermediate_df)
    
    # Look for mz or retention time columns to see if groups are ambiguous or not.
    # intermediate_df = pandas.DataFrame({'Metabolite':['name1', 'name2', 'name3', 'name4'], 'mz':[1,1, 2,2], 'rt':[0,0, 1,1]})
    # column_name_map = {column:column.lower() for column in intermediate_df.columns}
    # mz_column_names = column_finders['moverz_quant'].name_dict_match(column_name_map)
    # rt_column_names = column_finders['retention_time'].name_dict_match(column_name_map)
    # intermediate_df.loc[:, 'mz_matches'] = intermediate_df.apply(_mz_rt_match, args=(intermediate_df, mz_column_names), axis=1)
    # intermediate_df.loc[:, 'rt_matches'] = intermediate_df.apply(_mz_rt_match, args=(intermediate_df, rt_column_names), axis=1)
    # intermediate_df.loc[:, 'mz_rt_matches'] = intermediate_df.apply(lambda row: row['mz_matches'].intersection(row['rt_matches']), axis=1)
    # intermediate_df = _add_group_num(intermediate_df, 'mz_rt_group_num', 'mz_rt_matches')
    
    
    # intermediate_df.loc[:, 'duplicate_attributes'] = 'Unknown'
    # final_groups = intermediate_df.groupby('final_group_num', dropna=False)
    # for name, final_group in final_groups:
    #     has_iso_group = ~final_group['iso_group_num'].isna().all()
    #     has_fuzzy_group = ~final_group['fuzzy_group_num'].isna().all()
    #     has_mz_rt_group = ~final_group['fuzzy_group_num'].isna().all()
    #     mz_rt_NA = final_group['mz_rt_NA'].all()
        
    #     if (has_iso_group or has_fuzzy_group) and (mz_rt_NA or has_mz_rt_group):
    #         intermediate_df.loc[final_group.index, 'duplicate_attributes'] = 'Equivalent'
    #     if has_mz_rt_group and not (has_iso_group or has_fuzzy_group):
    #         intermediate_df.loc[final_group.index, 'duplicate_attributes'] = 'Ambiguous'
    
    # Create the DATA intermediate_df.
    data_intermediate_df = data_dup_group_df.merge(intermediate_df, on='Metabolite', how='left')
    names_in_data_not_met_mask = data_dup_group_df.loc[:, 'Metabolite'].isin(names_in_data_not_met)
    data_not_in_met_df = data_dup_group_df[names_in_data_not_met_mask]
    data_not_in_met_df = add_match_columns(data_not_in_met_df, column_finders)
    final_group_num = intermediate_df.loc[:, 'final_group_num'].max()+1
    final_group_num = final_group_num if not pandas.isna(final_group_num) else 0
    data_not_in_met_df = add_final_group_num(data_not_in_met_df, )
    data_not_in_met_df = add_duplicate_rank(data_not_in_met_df, ['str_len', 'Metabolite'])
    data_not_in_met_df = add_duplicate_attributes(data_not_in_met_df)
    data_not_in_met_df.loc[:, 'num_NA_values'] = pandas.NA
    data_intermediate_df = data_intermediate_df.loc[:, data_not_in_met_df.columns]
    data_intermediate_df.update(data_not_in_met_df)
    
    data_df_to_add = pandas.DataFrame({"duplicate_groups": [i]*data_intermediate_df.shape[0],
                                       "duplicate_sub_groups": data_intermediate_df.loc[:, 'final_group_num'].astype(int).values,
                                       "duplicate_ranks": data_intermediate_df.loc[:, 'Duplicate Rank'].values,
                                       "duplicate_attributes": data_intermediate_df.loc[:, 'duplicate_attributes'].astype(str).values},
                                      index = data_intermediate_df.index, dtype=str)
    data_group_dfs.append(data_df_to_add)
    
    met_df_to_add = pandas.DataFrame({"duplicate_groups": [i]*intermediate_df.shape[0],
                                      "duplicate_sub_groups": intermediate_df.loc[:, 'final_group_num'].values,
                                      "duplicate_ranks": intermediate_df.loc[:, 'Duplicate Rank'].values,
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
    
data_df.to_csv('C:/Users/Sparda/Desktop/Moseley Lab/Code/mwtab/scratch/temp1.csv', index=False)
met_df.to_csv('C:/Users/Sparda/Desktop/Moseley Lab/Code/mwtab/scratch/temp2.csv', index=False)


