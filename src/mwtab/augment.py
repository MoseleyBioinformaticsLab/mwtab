# -*- coding: utf-8 -*-
"""
Augment the DATA and METABOLITES tables of an mwtab file with extra columns.
"""

import pandas

from . import utility_functions

_create_numeric_df = utility_functions._create_numeric_df
find_metabolite_families = utility_functions.find_metabolite_families
_determine_row_subsets = utility_functions._determine_row_subsets


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




data_df = []
met_df = []


# Remove duplicates from Data.
data_duplicates = get_duplicate_rows(data_df, 'DATA')
met_duplicates = get_duplicate_rows(met_df, 'METABOLITES', False)


# Remove data_duplicates that are family members.
# data_metabolites = list(data_df.loc[:, "Metabolite"])
# met_to_root, root_to_mets = find_metabolite_families(data_metabolites)
# data_duplicate_mets = data_df.loc[data_duplicates.index, 'Metabolite']
# dup_indexes_to_keep = []
# if len(data_duplicates) != 0:
#     groups = data_duplicates.groupby(data_duplicates.columns.tolist(), dropna=False)
#     for name, group in groups:
#         names = list(data_df.loc[group.index, :].iloc[:, 0])
#         roots_to_names = {}
#         for name2 in names:
#             root = met_to_root[name2]
#             if root in roots_to_names:
#                 roots_to_names[root].append(name2)
#             else:
#                 roots_to_names[root] = [name2]
#         names_to_keep = [names[0] for root, names in roots_to_names.items() if len(names) == 1]
#         dup_indexes_to_keep += list(data_df[data_df.loc[:, 'Metabolite'].isin(names_to_keep)].index)
# data_duplicate_mets = data_duplicate_mets.loc[dup_indexes_to_keep]


# If there are duplicate measurements, we can remove from Data and 
# Metabolites if the row in Metabolites is all NA, or if the row 
# in Metabolites is duplicated.
met_to_root, root_to_mets = find_metabolite_families(list(data_df.loc[:, "Metabolite"]))
data_duplicate_mets = data_df.loc[data_duplicates.index, 'Metabolite']
dup_data_rows = met_df[met_df.loc[:, 'Metabolite'].isin(data_duplicate_mets)]
groups = data_duplicates.groupby(data_duplicates.columns.tolist(), dropna=False)
names_to_remove = []
names_in_data_not_met = []
for name, group in groups:
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
    
    temp_df = dup_data_rows[dup_data_rows.loc[:, 'Metabolite'].isin(names)]
    
    # If the names aren't in the METABOLITES table remove them from the DATA table.
    if len(temp_df.loc[:, 'Metabolite'].unique()) < len(set(names)):
        names_in_data_not_met += [name for name in names if name not in temp_df.loc[:, "Metabolite"].values]
        names = [name2 for name2 in names if name2 not in names_in_data_not_met]
    
    if len(temp_df) < 2:
        continue
    
    # temp_df = pandas.DataFrame([['asdf', pandas.NA, pandas.NA], 
    #                             ['qwer', 'a', 'b'], 
    #                             ['as', pandas.NA, pandas.NA ],
    #                             ['asd', 'a', pandas.NA],
    #                             ['qwe', 'a', 'b'],
    #                             ['zxcv', 'c', 'd']], columns=['Metabolite', 'col1', 'col2'])
    # names = list(temp_df.iloc[:,0])
    # names_to_remove = []
    
    # Find all NA rows in Metabolites.
    # all_na_rows_selector = temp_df.iloc[:, 1:].isna().all(axis=1)
    # if all_na_rows_selector.sum() == len(names):
    #     temp_dups_to_remove = temp_df[all_na_rows_selector].loc[:, ['Metabolite']]
    #     temp_dups_to_remove.loc[:, 'str_len'] = temp_dups_to_remove.loc[:, 'Metabolite'].str.len()
    #     temp_dups_to_remove = temp_dups_to_remove.sort_values(by=['str_len', 'Metabolite'], ascending=False)
    #     names_to_remove += list(temp_dups_to_remove.loc[:, 'Metabolite'].values[:-1])
    # else:
    #     names_to_remove += list(temp_df[all_na_rows_selector].loc[:, 'Metabolite'])
    # names = [name2 for name2 in names if name2 not in names_to_remove]
    # temp_df = temp_df[temp_df.loc[:, 'Metabolite'].isin(names)]
    # if len(temp_df) < 2:
    #     continue
    
    # Remove metabolites that have the same column info in the Metabolites table.
    met_group_duplicates_selector = met_duplicates.index.isin(temp_df.index)
    met_duplicate_indexes_to_remove = met_duplicates[met_group_duplicates_selector].index
    temp_dups_to_remove = temp_df.loc[met_duplicate_indexes_to_remove, ['Metabolite']]
    temp_dups_to_remove.loc[:, 'str_len'] = temp_dups_to_remove.loc[:, 'Metabolite'].str.len()
    temp_dups_to_remove = temp_dups_to_remove.sort_values(by=['str_len', 'Metabolite'])
    names_to_remove += list(temp_dups_to_remove.loc[:, 'Metabolite'].values[1:])
    names = [name2 for name2 in names if name2 not in names_to_remove]
    temp_df = temp_df[temp_df.loc[:, 'Metabolite'].isin(names)]
    if len(temp_df) < 2:
        continue
    
    # Remove row subsets from the Metabolites table.
    is_a_subset = _determine_row_subsets(temp_df.iloc[:, 1:])
    if sum(is_a_subset) == len(names):
        temp_dups_to_remove = temp_df[is_a_subset].loc[:, ['Metabolite']]
        temp_dups_to_remove.loc[:, 'str_len'] = temp_dups_to_remove.loc[:, 'Metabolite'].str.len()
        temp_dups_to_remove = temp_dups_to_remove.sort_values(by=['str_len', 'Metabolite'], ascending=False)
        names_to_remove += list(temp_dups_to_remove.loc[:, 'Metabolite'].values[:-1])
    else:
        names_to_remove += list(temp_df[is_a_subset].loc[:, 'Metabolite'])
    
if names_to_remove:
    met_df = met_df[~met_df.loc[:, 'Metabolite'].isin(names_to_remove)]
    tabfile.set_metabolites_from_pandas(met_df, data_section_key)
    data_df = data_df[~data_df.loc[:, 'Metabolite'].isin(names_to_remove)]
    tabfile.set_metabolites_data_from_pandas(data_df, data_section_key)
    print("The following metabolites in the 'DATA' and 'METABOLITES' sections were removed "
          "because they are duplicates of other metabolites:")
    names_to_remove = sorted(names_to_remove)
    try:
        print("\n".join([name for name in names_to_remove]))
    except UnicodeEncodeError:
        print("\n".join([name for name in names_to_remove]).encode('utf-8'))
if names_in_data_not_met:
    data_df = data_df[~data_df.loc[:, 'Metabolite'].isin(names_in_data_not_met)]
    tabfile.set_metabolites_data_from_pandas(data_df, data_section_key)
    print("The following metabolites in the 'DATA' section were removed "
          "because they are duplicates of other metabolites:")
    names_in_data_not_met = sorted(names_in_data_not_met)
    try:
        print("\n".join([name for name in names_in_data_not_met]))
    except UnicodeEncodeError:
        print("\n".join([name for name in names_in_data_not_met]).encode('utf-8'))




