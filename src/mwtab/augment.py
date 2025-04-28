# -*- coding: utf-8 -*-
"""
Augment the DATA and METABOLITES tables of an mwtab file with extra columns.
"""
import re

import pandas

from . import utility_functions

_create_numeric_df = utility_functions._create_numeric_df
find_metabolite_families = utility_functions.find_metabolite_families
_determine_row_subsets = utility_functions._determine_row_subsets
compute_fuzz_ratios = utility_functions.compute_fuzz_ratios
_compute_fuzz_ratio = utility_functions._compute_fuzz_ratio


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
    working_len = len(working_row['short_name'])
    for index, row in working_df.iterrows():
        iso_counts_match = working_row['iso_counts'] == row['iso_counts']
        name_in_name = working_row['short_name'] in row['short_name']
        names_fuzzy_match = row['fuzzy_matches'].get(working_row['Metabolite'])
        names_are_simliar_len = abs(working_len - len(row['short_name'])) < max_len_diff
        names_are_similar = ( names_are_simliar_len and name_in_name) or names_fuzzy_match
        if has_iso_count and iso_counts_match and names_are_similar:
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


GROUP_DEFAULT = 9999
def add_match_columns(df):
    """
    """
    df.loc[:, 'str_len'] = df.loc[:, 'Metabolite'].str.len()
    df.loc[:, 'iso_counts'] = df.loc[:, 'Metabolite'].apply(get_isotope_counts)
    df.loc[:, 'short_name'] = df.loc[:, 'Metabolite'].apply(remove_isotope_substring)
    df.loc[:, 'fuzzy_matches'] = df.apply(fuzzy_match, args=(df,), axis=1)
    df.loc[:, 'isotopologue_duplicates'] = df.apply(_isotopologue_match, args=(df, 4), axis=1)
    df.loc[:, 'isotopologue_duplicates'] = df.apply(_isotopologue_match_propagate, args=(df,), axis=1)
    name_to_iso_group = {}
    index_to_iso_group = {}
    name_to_fuzzy_group = {}
    index_to_fuzzy_group = {}
    iso_group_num = 0
    fuzzy_group_num = 0
    df.loc[:, 'iso_group_num'] = GROUP_DEFAULT
    df.loc[:, 'fuzzy_group_num'] = GROUP_DEFAULT
    for index, row in df.iterrows():
        name = row['Metabolite']
        if iso_dups := row['isotopologue_duplicates']:
            iso_group_found = False
            for iso_dup in iso_dups:
                if iso_dup in name_to_iso_group:
                    df.loc[index, 'iso_group_num'] = name_to_iso_group[iso_dup]
                    name_to_iso_group[name] = name_to_iso_group[iso_dup]
                    index_to_iso_group[index] = name_to_iso_group[iso_dup]
                    iso_group_found = True
            if not iso_group_found:
                df.loc[index, 'iso_group_num'] = iso_group_num
                name_to_iso_group[name] = iso_group_num
                index_to_iso_group[index] = iso_group_num
                for iso_dup in iso_dups:
                    name_to_iso_group[iso_dup] = iso_group_num
                    iso_dup_indexes = [iso_index for iso_index in df[df.loc[:, 'Metabolite'].isin([iso_dup])].index]
                    for iso_index in iso_dup_indexes:
                        index_to_iso_group[iso_index] = iso_group_num
                    df.loc[df[df.loc[:, 'Metabolite'].isin([iso_dup])].index, 'iso_group_num'] = iso_group_num
                iso_group_num += 1
        
        if fuzzy_dups := row['fuzzy_matches']:
            fuzzy_group_found = False
            for fuzzy_dup in fuzzy_dups:
                if fuzzy_dup in name_to_fuzzy_group:
                    df.loc[index, 'fuzzy_group_num'] = name_to_fuzzy_group[fuzzy_dup]
                    name_to_fuzzy_group[name] = name_to_fuzzy_group[fuzzy_dup]
                    index_to_fuzzy_group[index] = name_to_fuzzy_group[fuzzy_dup]
                    fuzzy_group_found = True
            if not fuzzy_group_found:
                df.loc[index, 'fuzzy_group_num'] = fuzzy_group_num
                name_to_fuzzy_group[name] = fuzzy_group_num
                index_to_fuzzy_group[index] = fuzzy_group_num
                for fuzzy_dup in fuzzy_dups:
                    name_to_fuzzy_group[fuzzy_dup] = fuzzy_group_num
                    fuzzy_dup_indexes = [fuzzy_index for fuzzy_index in df[df.loc[:, 'Metabolite'].isin([fuzzy_dup])].index]
                    for fuzzy_index in fuzzy_dup_indexes:
                        index_to_fuzzy_group[fuzzy_index] = fuzzy_group_num
                    df.loc[df[df.loc[:, 'Metabolite'].isin([fuzzy_dup])].index, 'fuzzy_group_num'] = fuzzy_group_num
                fuzzy_group_num += 1
    return df

def add_final_group_num(df, group_start_num=0):
    """
    """
    final_group_num = group_start_num
    df.loc[:, 'final_group_num'] = GROUP_DEFAULT
    iso_group_nums = {num for num in df.loc[:, 'iso_group_num'].values if num != GROUP_DEFAULT}
    iso_group_nums = sorted(list(iso_group_nums))
    for iso_group_num in iso_group_nums:
        iso_group_mask = df.loc[:, 'iso_group_num'].isin([iso_group_num])
        final_group_mask = df.loc[:, 'final_group_num'].isin([GROUP_DEFAULT])
        group_df = df[(iso_group_mask & final_group_mask)]
        info_group_nums = {num for info_set in group_df.loc[:, 'info_match_group'].values for num in info_set}
        info_group_nums = sorted(list(info_group_nums))
        for info_group_num in info_group_nums:
            info_group_mask = group_df.loc[:, 'info_match_group'].apply(lambda x: info_group_num in x)
            sub_final_group_mask = df.loc[group_df.index, 'final_group_num'].isin([GROUP_DEFAULT])
            final_group_df = group_df[(info_group_mask & sub_final_group_mask)]
            if final_group_df.shape[0]:
                df.loc[final_group_df.index, 'final_group_num'] = final_group_num
                final_group_num += 1
    
    fuzzy_group_nums = {num for num in df.loc[:, 'fuzzy_group_num'].values if num != GROUP_DEFAULT}
    fuzzy_group_nums = sorted(list(fuzzy_group_nums))
    for fuzzy_group_num in fuzzy_group_nums:
        fuzzy_group_mask = df.loc[:, 'fuzzy_group_num'].isin([fuzzy_group_num])
        final_group_mask = df.loc[:, 'final_group_num'].isin([GROUP_DEFAULT])
        group_df = df[(fuzzy_group_mask & final_group_mask)]
        info_group_nums = {num for info_set in group_df.loc[:, 'info_match_group'].values for num in info_set}
        info_group_nums = sorted(list(info_group_nums))
        for info_group_num in info_group_nums:
            info_group_mask = group_df.loc[:, 'info_match_group'].apply(lambda x: info_group_num in x)
            sub_final_group_mask = df.loc[group_df.index, 'final_group_num'].isin([GROUP_DEFAULT])
            final_group_df = group_df[(info_group_mask & sub_final_group_mask)]
            if final_group_df.shape[0]:
                df.loc[final_group_df.index, 'final_group_num'] = final_group_num
                final_group_num += 1
    
    # Set the rest to their own groups.
    final_group_mask = df.loc[:, 'final_group_num'].isin([GROUP_DEFAULT])
    for index, row in df[final_group_mask].iterrows():
        df.loc[index, 'final_group_num'] = final_group_num
        final_group_num += 1
    
    return df

def set_duplicate_rank(df, sort_by):
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


"""
I am thinking detecting duplicates should add 3 columns. 
1 like "Duplicate_Groups" that assigns a unique number to each group of duplicates.
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

met_df = pandas.DataFrame({'Metabolite': ['alanine-13C3-15N', 'L-ALANINE_13C3', 'alanine-13c3-15n0', 'r-alanine', 'r_alanine', 'alanine', 'not-alanine', 'NOT-alanine', 'also_not_alanine'], 
                           'col1': ['1', '1', pandas.NA, '1', pandas.NA, pandas.NA, '2', '2', '1'],
                           'col2': ['2', '2', pandas.NA, pandas.NA, pandas.NA, pandas.NA, '3', '3', '4']})

data_df = pandas.DataFrame({'Metabolite': ['alanine-13C3-15N', 'L-ALANINE_13C3', 'alanine-13c3-15n0', 'r-alanine', 'r_alanine', 'alanine', 'not-alanine', 'NOT-alanine', 'not_in_met', 'also_not_alanine'], 
                           'col1': ['1', '1', '1', '1', '1', '1', '1', '1', '1', '1'],
                           'col2': ['2', '2', '2', '2', '2', '2', '2', '2', '2', '2']})


# Remove duplicates from Data.
data_duplicates = get_duplicate_rows(data_df, 'DATA')
met_duplicates = get_duplicate_rows(met_df, 'METABOLITES', False)


# If there are duplicate measurements, we can remove from Data and 
# Metabolites if the row in Metabolites is all NA, or if the row 
# in Metabolites is duplicated.
met_to_root, root_to_mets = find_metabolite_families(list(data_df.loc[:, "Metabolite"]))
data_duplicate_mets = data_df.loc[data_duplicates.index, 'Metabolite']
met_rows_of_data_dups = met_df[met_df.loc[:, 'Metabolite'].isin(data_duplicate_mets)]
groups = data_duplicates.groupby(data_duplicates.columns.tolist(), dropna=False)
names_to_remove = []
data_group_indexes = []
met_group_indexes = []
data_dup_rank_indexes = []
met_dup_rank_indexes = []
# These are category 1, duplicates in DATA where the metabolite is not in METABOLITES. Very strong removal recomendation.
names_in_data_not_met = []
indexes_in_data_not_met = []
# Category 2, duplicates that have the same information in each column. Recommended removal, but which one to remove needs to be determined.
names_in_data_with_identical_info_in_met = []
indexes_in_data_with_identical_info_in_met = []
names_in_met_with_identical_info_in_met = []
indexes_in_met_with_identical_info_in_met = []
# Category 3, duplicates that have a subset of the information in each column. Recommended removal.
names_in_data_with_subset_info_in_met = []
indexes_in_data_with_subset_info_in_met = []
names_in_met_with_subset_info_in_met = []
indexes_in_met_with_subset_info_in_met = []
# Store data frames for each gorup to update the main dataframes with.
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
    data_group_indexes.append(data_dup_group_df.index)
    
    met_dup_group_df = met_rows_of_data_dups[met_rows_of_data_dups.loc[:, 'Metabolite'].isin(names)]
    met_group_indexes.append(met_dup_group_df.index)
    intermediate_df = met_dup_group_df.copy()
    
    # If the names aren't in the METABOLITES table remove them from the DATA table.
    if len(met_dup_group_df.loc[:, 'Metabolite'].unique()) < len(set(names)):
        # Note, if somehow all the names aren't in METABOLITES, this will add them all to the list.
        names_in_data_not_met += [name for name in names if name not in met_dup_group_df.loc[:, "Metabolite"].values]
        # names = [name2 for name2 in names if name2 not in names_in_data_not_met]
        indexes_in_data_not_met += list(data_dup_group_df.loc[:, 'Metabolite'].isin(names_in_data_not_met).index)
        # data_dup_group_df = data_dup_group_df[data_dup_group_df.loc[:, 'Metabolite'].isin(names)]
    
    
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
    
    intermediate_df = add_match_columns(intermediate_df)
    intermediate_df.loc[:, 'num_NA_values'] = intermediate_df.isna().sum(axis=1)
    
    # intermediate_df.loc[:, 'str_len'] = intermediate_df.loc[:, 'Metabolite'].str.len()
    # intermediate_df.loc[:, 'num_NA_values'] = intermediate_df.isna().sum(axis=1)
    # intermediate_df.loc[:, 'iso_counts'] = intermediate_df.loc[:, 'Metabolite'].apply(get_isotope_counts)
    # intermediate_df.loc[:, 'short_name'] = intermediate_df.loc[:, 'Metabolite'].apply(remove_isotope_substring)
    # intermediate_df.loc[:, 'short_full_name'] = intermediate_df.loc[:, 'Metabolite'].str.lower()
    # intermediate_df.loc[:, 'fuzzy_matches'] = intermediate_df.apply(fuzzy_match, args=(intermediate_df,), axis=1)
    # intermediate_df.loc[:, 'isotopologue_duplicates'] = intermediate_df.apply(_isotopologue_match, args=(intermediate_df, 4), axis=1)
    # intermediate_df.loc[:, 'isotopologue_duplicates'] = intermediate_df.apply(_isotopologue_match_propagate, args=(intermediate_df,), axis=1)
    # isotopologue_names = [row['Metabolite'] for index, row in intermediate_df.iterrows() if row['isotopologue_duplicates']]
    # name_to_iso_group = {}
    # iso_group_to_name = {}
    # index_to_iso_group = {}
    # name_to_fuzzy_group = {}
    # index_to_fuzzy_group = {}
    # iso_group_num = 0
    # fuzzy_group_num = 0
    # GROUP_DEFAULT = 9999
    # intermediate_df.loc[:, 'iso_group_num'] = GROUP_DEFAULT
    # intermediate_df.loc[:, 'fuzzy_group_num'] = GROUP_DEFAULT
    # for index, row in intermediate_df.iterrows():
    #     name = row['Metabolite']
    #     if iso_dups := row['isotopologue_duplicates']:
    #         iso_group_found = False
    #         for iso_dup in iso_dups:
    #             if iso_dup in name_to_iso_group:
    #                 intermediate_df.loc[index, 'iso_group_num'] = name_to_iso_group[iso_dup]
    #                 name_to_iso_group[name] = name_to_iso_group[iso_dup]
    #                 index_to_iso_group[index] = name_to_iso_group[iso_dup]
    #                 iso_group_found = True
    #         if not iso_group_found:
    #             intermediate_df.loc[index, 'iso_group_num'] = iso_group_num
    #             name_to_iso_group[name] = iso_group_num
    #             index_to_iso_group[index] = iso_group_num
    #             for iso_dup in iso_dups:
    #                 name_to_iso_group[iso_dup] = iso_group_num
    #                 iso_dup_indexes = [iso_index for iso_index in intermediate_df[intermediate_df.loc[:, 'Metabolite'].isin([iso_dup])].index]
    #                 for iso_index in iso_dup_indexes:
    #                     index_to_iso_group[iso_index] = iso_group_num
    #                 intermediate_df.loc[intermediate_df[intermediate_df.loc[:, 'Metabolite'].isin([iso_dup])].index, 'iso_group_num'] = iso_group_num
    #             iso_group_num += 1
        
    #     if fuzzy_dups := row['fuzzy_matches']:
    #         fuzzy_group_found = False
    #         for fuzzy_dup in fuzzy_dups:
    #             if fuzzy_dup in name_to_fuzzy_group:
    #                 intermediate_df.loc[index, 'fuzzy_group_num'] = name_to_fuzzy_group[fuzzy_dup]
    #                 name_to_fuzzy_group[name] = name_to_fuzzy_group[fuzzy_dup]
    #                 index_to_fuzzy_group[index] = name_to_fuzzy_group[fuzzy_dup]
    #                 fuzzy_group_found = True
    #         if not fuzzy_group_found:
    #             intermediate_df.loc[index, 'fuzzy_group_num'] = fuzzy_group_num
    #             name_to_fuzzy_group[name] = fuzzy_group_num
    #             index_to_fuzzy_group[index] = fuzzy_group_num
    #             for fuzzy_dup in fuzzy_dups:
    #                 name_to_fuzzy_group[fuzzy_dup] = fuzzy_group_num
    #                 fuzzy_dup_indexes = [fuzzy_index for fuzzy_index in intermediate_df[intermediate_df.loc[:, 'Metabolite'].isin([fuzzy_dup])].index]
    #                 for fuzzy_index in fuzzy_dup_indexes:
    #                     index_to_fuzzy_group[fuzzy_index] = fuzzy_group_num
    #                 intermediate_df.loc[intermediate_df[intermediate_df.loc[:, 'Metabolite'].isin([fuzzy_dup])].index, 'fuzzy_group_num'] = fuzzy_group_num
    #             fuzzy_group_num += 1
    
    intermediate_df = add_final_group_num(intermediate_df, group_start_num=0)
    
    # intermediate_df.loc[:, 'final_group_num'] = GROUP_DEFAULT
    # iso_group_nums = {num for num in intermediate_df.loc[:, 'iso_group_num'].values if num != GROUP_DEFAULT}
    # iso_group_nums = sorted(list(iso_group_nums))
    # final_group_num = 0
    # for iso_group_num in iso_group_nums:
    #     iso_group_mask = intermediate_df.loc[:, 'iso_group_num'].isin([iso_group_num])
    #     final_group_mask = intermediate_df.loc[:, 'final_group_num'].isin([GROUP_DEFAULT])
    #     group_df = intermediate_df[(iso_group_mask & final_group_mask)]
    #     info_group_nums = {num for info_set in group_df.loc[:, 'info_match_group'].values for num in info_set}
    #     info_group_nums = sorted(list(info_group_nums))
    #     for info_group_num in info_group_nums:
    #         info_group_mask = group_df.loc[:, 'info_match_group'].apply(lambda x: info_group_num in x)
    #         sub_final_group_mask = intermediate_df.loc[group_df.index, 'final_group_num'].isin([GROUP_DEFAULT])
    #         final_group_df = group_df[(info_group_mask & sub_final_group_mask)]
    #         if final_group_df.shape[0]:
    #             intermediate_df.loc[final_group_df.index, 'final_group_num'] = final_group_num
    #             final_group_num += 1
    
    # fuzzy_group_nums = {num for num in intermediate_df.loc[:, 'fuzzy_group_num'].values if num != GROUP_DEFAULT}
    # fuzzy_group_nums = sorted(list(fuzzy_group_nums))
    # for fuzzy_group_num in fuzzy_group_nums:
    #     fuzzy_group_mask = intermediate_df.loc[:, 'fuzzy_group_num'].isin([fuzzy_group_num])
    #     final_group_mask = intermediate_df.loc[:, 'final_group_num'].isin([GROUP_DEFAULT])
    #     group_df = intermediate_df[(fuzzy_group_mask & final_group_mask)]
    #     info_group_nums = {num for info_set in group_df.loc[:, 'info_match_group'].values for num in info_set}
    #     info_group_nums = sorted(list(info_group_nums))
    #     for info_group_num in info_group_nums:
    #         info_group_mask = group_df.loc[:, 'info_match_group'].apply(lambda x: info_group_num in x)
    #         sub_final_group_mask = intermediate_df.loc[group_df.index, 'final_group_num'].isin([GROUP_DEFAULT])
    #         final_group_df = group_df[(info_group_mask & sub_final_group_mask)]
    #         if final_group_df.shape[0]:
    #             intermediate_df.loc[final_group_df.index, 'final_group_num'] = final_group_num
    #             final_group_num += 1
    
    # # Set the rest to their own groups.
    # final_group_mask = intermediate_df.loc[:, 'final_group_num'].isin([GROUP_DEFAULT])
    # for index, row in intermediate_df[final_group_mask].iterrows():
    #     intermediate_df.loc[index, 'final_group_num'] = final_group_num
    #     final_group_num += 1
    
    
    intermediate_df = set_duplicate_rank(intermediate_df, ['num_NA_values', 'str_len', 'Metabolite'])
    # Set the rank for the duplicates in each group.
    # intermediate_df.loc[:, 'Duplicate Rank'] = ''
    # final_group_nums = {num for num in intermediate_df.loc[:, 'final_group_num'].values}
    # final_group_nums = sorted(list(final_group_nums))
    # for final_group_num in final_group_nums:
    #     group_df = intermediate_df[intermediate_df.loc[:, 'final_group_num'].isin([final_group_num])]
    #     group_df = group_df.sort_values(by=['num_NA_values', 'str_len', 'Metabolite'])
    #     group_df.loc[:, 'Duplicate Rank'] = [str(num) for num in range(group_df.shape[0])]
    #     intermediate_df.update(group_df.loc[:, 'Duplicate Rank'])
    
    # Create the DATA intermediate_df.
    data_intermediate_df = data_dup_group_df.merge(intermediate_df, on='Metabolite', how='left')
    names_in_data_not_met_mask = data_dup_group_df.loc[:, 'Metabolite'].isin(names_in_data_not_met)
    data_not_in_met_df = data_dup_group_df[names_in_data_not_met_mask]
    data_not_in_met_df = add_match_columns(data_not_in_met_df)
    data_not_in_met_df = add_final_group_num(data_not_in_met_df, intermediate_df.loc[:, 'final_group_num'].max()+1)
    data_not_in_met_df = set_duplicate_rank(data_not_in_met_df, ['str_len', 'Metabolite'])
    data_not_in_met_df.loc[:, 'num_NA_values'] = pandas.NA
    data_intermediate_df = data_intermediate_df.loc[:, data_not_in_met_df.columns]
    data_intermediate_df.update(data_not_in_met_df)
    
    data_df_to_add = pandas.DataFrame({"Duplicate_Groups": [i]*data_intermediate_df.shape[0],
                                       "Duplicate_Sub_Groups": data_intermediate_df.loc[:, 'final_group_num'].astype(int).values,
                                       "Duplicate_Ranks": data_intermediate_df.loc[:, 'Duplicate Rank'].values},
                                      index = data_intermediate_df.index, dtype=str)
    data_group_dfs.append(data_df_to_add)
    
    met_df_to_add = pandas.DataFrame({"Duplicate_Groups": [i]*intermediate_df.shape[0],
                                      "Duplicate_Sub_Groups": intermediate_df.loc[:, 'final_group_num'].values,
                                      "Duplicate_Ranks": intermediate_df.loc[:, 'Duplicate Rank'].values},
                                     index = intermediate_df.index, dtype=str)
    met_group_dfs.append(met_df_to_add)
        
    
    
    
    # if len(met_dup_group_df) > 1:
    #     # Find metabolites that have the same column info in the Metabolites table.
    #     met_group_duplicates_selector = met_duplicates.index.isin(met_dup_group_df.index)
    #     met_duplicate_indexes_to_remove = met_duplicates[met_group_duplicates_selector].index
    #     indexes_in_met_with_identical_info_in_met += list(met_duplicate_indexes_to_remove)
    #     met_names = met_dup_group_df.loc[met_duplicate_indexes_to_remove, 'Metabolite']
    #     names_in_met_with_identical_info_in_met += list(met_names)
    #     data_indexes = list(data_dup_group_df[data_dup_group_df.loc[:, 'Metabolite'].isin(names_in_met_with_identical_info_in_met)].index)
    #     indexes_in_data_with_identical_info_in_met += data_indexes
    #     names_in_data_with_identical_info_in_met += list(data_dup_group_df.loc[data_indexes, 'Metabolite'])
        
        
    #     # Look to see if the names isotopologue counts match and if the names fuzzy match strongly.
    #     # If so, we can classify them as likely direct duplicates.
    #     iso_df = met_dup_group_df.copy()
    #     iso_df.loc[:, 'str_len'] = iso_df.loc[:, 'Metabolite'].str.len()
    #     iso_df.loc[:, 'num_NA_values'] = iso_df.isna().sum(axis=1)
    #     iso_df.loc[:, 'iso_counts'] = iso_df.loc[:, 'Metabolite'].apply(get_isotope_counts)
    #     iso_df.loc[:, 'short_name'] = iso_df.loc[:, 'Metabolite'].apply(remove_isotope_substring)
    #     iso_df.loc[:, 'short_full_name'] = iso_df.loc[:, 'Metabolite'].str.lower()
    #     short_to_long = {record['short_name']:record['Metabolite'] for record in iso_df.to_dict('records')}
    #     iso_df.loc[:, 'fuzzy_matches'] = iso_df.apply(fuzzy_match, args=(iso_df,), axis=1)
    #     # iso_df.loc[:, 'fuzzy_matches'] = iso_df.loc[:, 'short_name'].apply(fuzzy_match, args=(iso_df.loc[:, 'short_name'].values, short_to_long))
    #     iso_df.loc[:, 'isotopologue_duplicates'] = iso_df.apply(_isotopologue_match, args=(iso_df, 4), axis=1)
    #     iso_df.loc[:, 'isotopologue_duplicates'] = iso_df.apply(_isotopologue_match_propagate, args=(iso_df,), axis=1)
    #     isotopologue_names = [row['Metabolite'] for index, row in iso_df.iterrows() if row['isotopologue_duplicates']]
    #     name_to_iso_group = {}
    #     iso_group_to_name = {}
    #     index_to_iso_group = {}
    #     name_to_fuzzy_group = {}
    #     index_to_fuzzy_group = {}
    #     iso_group_num = 0
    #     fuzzy_group_num = 0
    #     GROUP_DEFAULT = 9999
    #     iso_df.loc[:, 'iso_group_num'] = GROUP_DEFAULT
    #     iso_df.loc[:, 'fuzzy_group_num'] = GROUP_DEFAULT
    #     for index, row in iso_df.iterrows():
    #         name = row['Metabolite']
    #         if iso_dups := row['isotopologue_duplicates']:
    #             iso_group_found = False
    #             for iso_dup in iso_dups:
    #                 if iso_dup in name_to_iso_group:
    #                     iso_df.loc[index, 'iso_group_num'] = name_to_iso_group[iso_dup]
    #                     name_to_iso_group[name] = name_to_iso_group[iso_dup]
    #                     index_to_iso_group[index] = name_to_iso_group[iso_dup]
    #                     iso_group_found = True
    #             if not iso_group_found:
    #                 iso_df.loc[index, 'iso_group_num'] = iso_group_num
    #                 name_to_iso_group[name] = iso_group_num
    #                 index_to_iso_group[index] = iso_group_num
    #                 for iso_dup in iso_dups:
    #                     name_to_iso_group[iso_dup] = iso_group_num
    #                     iso_dup_indexes = [iso_index for iso_index in iso_df[iso_df.loc[:, 'Metabolite'].isin([iso_dup])].index]
    #                     for iso_index in iso_dup_indexes:
    #                         index_to_iso_group[iso_index] = iso_group_num
    #                     iso_df.loc[iso_df[iso_df.loc[:, 'Metabolite'].isin([iso_dup])].index, 'iso_group_num'] = iso_group_num
    #                 iso_group_num += 1
            
    #         if fuzzy_dups := row['fuzzy_matches']:
    #             fuzzy_group_found = False
    #             for fuzzy_dup in fuzzy_dups:
    #                 if fuzzy_dup in name_to_fuzzy_group:
    #                     iso_df.loc[index, 'fuzzy_group_num'] = name_to_fuzzy_group[fuzzy_dup]
    #                     name_to_fuzzy_group[name] = name_to_fuzzy_group[fuzzy_dup]
    #                     index_to_fuzzy_group[index] = name_to_fuzzy_group[fuzzy_dup]
    #                     fuzzy_group_found = True
    #             if not fuzzy_group_found:
    #                 iso_df.loc[index, 'fuzzy_group_num'] = fuzzy_group_num
    #                 name_to_fuzzy_group[name] = fuzzy_group_num
    #                 index_to_fuzzy_group[index] = fuzzy_group_num
    #                 for fuzzy_dup in fuzzy_dups:
    #                     name_to_fuzzy_group[fuzzy_dup] = fuzzy_group_num
    #                     fuzzy_dup_indexes = [fuzzy_index for fuzzy_index in iso_df[iso_df.loc[:, 'Metabolite'].isin([fuzzy_dup])].index]
    #                     for fuzzy_index in fuzzy_dup_indexes:
    #                         index_to_fuzzy_group[fuzzy_index] = fuzzy_group_num
    #                     iso_df.loc[iso_df[iso_df.loc[:, 'Metabolite'].isin([fuzzy_dup])].index, 'fuzzy_group_num'] = fuzzy_group_num
    #                 fuzzy_group_num += 1
        
    #     #TODO Create a duplicates ranking within iso groups and within fuzzy groups.
    #     iso_groups = iso_df.groupby('iso_group_num')
    #     for iso_name, iso_group in iso_groups:
    #         print(iso_name)
        
    #     # met_duplicate_indexes_to_remove = met_duplicates[met_group_duplicates_selector].index
    #     # temp_dups_to_remove = met_dup_group_df.loc[met_duplicate_indexes_to_remove, ['Metabolite']]
    #     # temp_dups_to_remove.loc[:, 'str_len'] = temp_dups_to_remove.loc[:, 'Metabolite'].str.len()
    #     # temp_dups_to_remove = temp_dups_to_remove.sort_values(by=['str_len', 'Metabolite'])
    #     # names_to_remove += list(temp_dups_to_remove.loc[:, 'Metabolite'].values[1:])
        
    #     # names = [name2 for name2 in names if name2 not in names_in_data_with_identical_info_in_met]
    #     # met_dup_group_df = met_dup_group_df[met_dup_group_df.loc[:, 'Metabolite'].isin(names)]
    #     # data_dup_group_df = data_dup_group_df[data_dup_group_df.loc[:, 'Metabolite'].isin(names)]
    #     if len(met_dup_group_df) > 1:
    #         # Find metabolites that have a subset of column information in the Metabolites table.
    #         is_a_subset = _determine_row_subsets(met_dup_group_df.iloc[:, 1:])
    #         # if sum(is_a_subset) == len(names):
    #         #     temp_dups_to_remove = met_dup_group_df[is_a_subset].loc[:, ['Metabolite']]
    #         #     temp_dups_to_remove.loc[:, 'str_len'] = temp_dups_to_remove.loc[:, 'Metabolite'].str.len()
    #         #     temp_dups_to_remove = temp_dups_to_remove.sort_values(by=['str_len', 'Metabolite'], ascending=False)
    #         #     names_to_remove += list(temp_dups_to_remove.loc[:, 'Metabolite'].values[:-1])
    #         # else:
    #         #     names_to_remove += list(met_dup_group_df[is_a_subset].loc[:, 'Metabolite'])
                
    #         names_in_met_with_subset_info_in_met += list(met_dup_group_df[is_a_subset].loc[:, 'Metabolite'])
    #         indexes_in_met_with_subset_info_in_met += list(met_dup_group_df[is_a_subset].index)
    #         data_indexes = list(data_dup_group_df[data_dup_group_df.loc[:, 'Metabolite'].isin(names_in_met_with_subset_info_in_met)].index)
    #         indexes_in_data_with_subset_info_in_met += data_indexes
    #         names_in_data_with_subset_info_in_met += list(data_dup_group_df.loc[data_indexes, 'Metabolite'])
            
    # met_dup_group_df = met_df.loc[met_group_indexes[-1], :]
    # data_dup_group_df = data_df.loc[data_group_indexes[-1], ['Metabolite']]
    # # met_dup_group_df = met_rows_of_data_dups[met_rows_of_data_dups.loc[:, 'Metabolite'].isin(original_names)]
    # # data_dup_group_df = data_df.loc[group.index, :].loc[:, ['Metabolite']]
    # # data_dup_group_df = data_dup_group_df[data_dup_group_df.loc[:, 'Metabolite'].isin(original_names)]
    
    # met_dup_group_df.loc[:, 'num_NA_values'] = met_dup_group_df.isna().sum(axis=1)
    # met_dup_group_df.loc[:, 'is_not_iso_dup'] = ~met_dup_group_df.loc[:, 'Metabolite'].isin(isotopologue_names)
    # met_dup_group_df.loc[:, 'has_identical_info'] = met_dup_group_df.loc[:, 'Metabolite'].isin(names_in_met_with_identical_info_in_met)
    # met_dup_group_df.loc[:, 'has_subset_info'] = met_dup_group_df.loc[:, 'Metabolite'].isin(names_in_met_with_subset_info_in_met)
    # met_dup_group_df.loc[:, 'str_len'] = met_dup_group_df.loc[:, 'Metabolite'].str.len()
    # met_dup_group_df = met_dup_group_df.sort_values(by=['num_NA_values', 'is_not_iso_dup', 'has_identical_info', 'has_subset_info', 'str_len', 'Metabolite'])
    # # met_dup_group_df.loc[:, 'Duplicate_Rank'] = range(1, met_dup_group_df.shape[0]+1)
    # met_dup_rank_indexes.append(met_dup_group_df.index)
    # met_dup_names = [row['Metabolite'] for index, row in met_dup_group_df.iterrows()]
    
    # data_not_in_met_df = data_dup_group_df[data_dup_group_df.loc[:, 'Metabolite'].isin(names_in_data_not_met)]
    # data_not_in_met_df.loc[:, 'str_len'] = data_not_in_met_df.loc[:, 'Metabolite'].str.len()
    # data_not_in_met_df = data_not_in_met_df.sort_values(by=['str_len', 'Metabolite'])
    # # data_dup_group_df.loc[:, 'Duplicate_Rank'] = range(1, data_dup_group_df.shape[0]+1)
    # data_indexes = [index for met_name in met_dup_names for index in data_dup_group_df[data_dup_group_df.loc[:, 'Metabolite'].isin([met_name])].index]
    # data_indexes += [row.name for index, row in data_not_in_met_df.iterrows()]
    # data_dup_rank_indexes.append(data_indexes)


if data_group_indexes:
    # TODO if Sub_Groups are identical to Ranks or Duplicate_Groups, then drop the column.
    data_df.loc[:, 'Duplicate_Groups'] = ''
    data_df.loc[:, 'Duplicate_Sub_Groups'] = ''
    data_df.loc[:, 'Duplicate_Ranks'] = ''
    for df in data_group_dfs:
        data_df.update(df)
    
    
    
    # data_df.loc[:, 'Duplicate_Groups'] = ''
    # for i, indexes in enumerate(data_group_indexes):
    #     data_df.loc[indexes, 'Duplicate_Groups'] = str(i)
    # data_df.loc[:, 'Duplicate_Rank'] = ''
    # for indexes in data_dup_rank_indexes:
    #     for i, index in enumerate(indexes):
    #         data_df.loc[index, 'Duplicate_Rank'] = str(i+1)
if met_group_indexes:
    met_df.loc[:, 'Duplicate_Groups'] = ''
    met_df.loc[:, 'Duplicate_Sub_Groups'] = ''
    met_df.loc[:, 'Duplicate_Ranks'] = ''
    for df in met_group_dfs:
        met_df.update(df)
    
    
    # met_df.loc[:, 'Duplicate_Groups'] = ''
    # for i, indexes in enumerate(met_group_indexes):
    #     met_df.loc[indexes, 'Duplicate_Groups'] = str(i)
    # met_df.loc[:, 'Duplicate_Rank'] = ''
    # for indexes in met_dup_rank_indexes:
    #     for i, index in enumerate(indexes):
    #         met_df.loc[index, 'Duplicate_Rank'] = str(i+1)

if indexes_in_data_not_met or indexes_in_data_with_identical_info_in_met or indexes_in_data_with_subset_info_in_met:
    data_df.loc[:, 'Duplicate_Info'] = ''
    data_df.loc[indexes_in_data_not_met, 'Duplicate_Info'] = 'Metabolite not in the METABOLITES section/table.'
    message = 'Metabolite metadata in the METABOLITES section/table has identical metadata with at least one of its duplicates.'
    data_df.loc[indexes_in_data_with_identical_info_in_met, 'Duplicate_Info'] = message
    message = 'Metabolite metadata in the METABOLITES section/table is a subset of the metadata for at least one of its duplicates.'
    data_df.loc[indexes_in_data_with_subset_info_in_met, 'Duplicate_Info'] = message

if indexes_in_met_with_identical_info_in_met or indexes_in_met_with_subset_info_in_met:
    met_df.loc[:, 'Duplicate_Info'] = ''
    message = 'Metabolite metadata in the METABOLITES section/table has identical metadata with at least one of its duplicates.'
    met_df.loc[indexes_in_met_with_identical_info_in_met, 'Duplicate_Info'] = message
    message = 'Metabolite metadata in the METABOLITES section/table is a subset of the metadata for at least one of its duplicates.'
    met_df.loc[indexes_in_met_with_subset_info_in_met, 'Duplicate_Info'] = message
    
# if names_to_remove:
#     met_df = met_df[~met_df.loc[:, 'Metabolite'].isin(names_to_remove)]
#     tabfile.set_metabolites_from_pandas(met_df, data_section_key)
#     data_df = data_df[~data_df.loc[:, 'Metabolite'].isin(names_to_remove)]
#     tabfile.set_metabolites_data_from_pandas(data_df, data_section_key)
#     print("The following metabolites in the 'DATA' and 'METABOLITES' sections were removed "
#           "because they are duplicates of other metabolites:")
#     names_to_remove = sorted(names_to_remove)
#     try:
#         print("\n".join([name for name in names_to_remove]))
#     except UnicodeEncodeError:
#         print("\n".join([name for name in names_to_remove]).encode('utf-8'))
# if names_in_data_not_met:
#     data_df = data_df[~data_df.loc[:, 'Metabolite'].isin(names_in_data_not_met)]
#     tabfile.set_metabolites_data_from_pandas(data_df, data_section_key)
#     print("The following metabolites in the 'DATA' section were removed "
#           "because they are duplicates of other metabolites:")
#     names_in_data_not_met = sorted(names_in_data_not_met)
#     try:
#         print("\n".join([name for name in names_in_data_not_met]))
#     except UnicodeEncodeError:
#         print("\n".join([name for name in names_in_data_not_met]).encode('utf-8'))




