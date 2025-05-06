# -*- coding: utf-8 -*-
"""
Code that repairs shifted rows in METABOLITES tables.
"""

import re
from typing import Tuple

import pandas



WRAP_STRING = r'[^a-zA-Z0-9]'
def _create_column_regex_any(match_strings: list[str]) -> str:
    """Return a regular expression string that will match any of the strings in match_strings.
    
    Automatically adds a WRAP_STRING to either side of each string to help with some fuzzy matching. 
    Intended to be used to match column names in METABOLITES data.
    
    Args:
        match_strings: list of strings to create a regex for.
    
    Returns:
        A regular expression string that will match any of the strings in match_strings.
    """
    regex = '|'.join(['(' + WRAP_STRING + '|^)' + match_string + '(' + WRAP_STRING + '|$)' for match_string in match_strings])
    return regex

def _create_column_regex_all(match_strings_sets: list[list[str]]) -> str:
    """Return a regular expression string that will match any of the string sets in match_strings_sets.
    
    Automatically adds a WRAP_STRING to either side of each string to help with some fuzzy matching. 
    Intended to be used to match column names in METABOLITES data. The "all" refers to requiring each 
    string within a set to be present to match, but if multiple sets are given, then any set will match. 
    "set" does not mean the actual "set" data type, any iterable collection is fine.
    
    Args:
        match_strings_sets: list of string sets to create a regex for.
    
    Returns:
        A regular expression string that will match any of the string sets in match_strings_sets.
    """
    regex = '|'.join([''.join(['(?=.*(' + WRAP_STRING + '|^)' + match_string + '(' + WRAP_STRING + '|$))' for match_string in match_strings]) for match_strings in match_strings_sets])
    return regex


class NameMatcher():
    def __init__(self, regex_search_strings: None|list[str] = None, 
                 not_regex_search_strings: None|list[str] = None, regex_search_sets: None|list[list[str]] = None,
                 in_strings: None|list[str] = None, not_in_strings: None|list[str] = None, 
                 in_string_sets: None|list[list[str]] = None, exact_strings: None|list[str] = None):
        self.regex_search_strings = regex_search_strings if regex_search_strings else []
        self.not_regex_search_strings = not_regex_search_strings if not_regex_search_strings else []
        self.regex_search_sets = regex_search_sets if regex_search_sets else []
        self.in_strings = in_strings if in_strings else []
        self.not_in_strings = not_in_strings if not_in_strings else []
        self.in_string_sets = in_string_sets if in_string_sets else []
        self.exact_strings = exact_strings if exact_strings else []
        
    def dict_match(self, name_map: dict[str, str]) -> list[str]:
        """Return a list of names that match based on the NameMatcher attributes.
        
        Find all names in name_map that match. name_map should be a dictionary of original names 
        to modified names. The value is used for matching, but the key is what will be returned. 
        Each of the name regex, in_string, and exact strings attributes are ORed together, 
        meaning any of them can be used to match, except for the "not" parameters. If a column 
        name is matched by a "not" parameter, then it overrides other matches and will be filtered out.
        
        Args:
            name_map: a dictionary of original names to the modified version of that name to use for matching.
        
        Returns:
            A list of names that match based on the NameMatcher attributes.
        """
        search_regex = _create_column_regex_any(self.regex_search_strings)
        search_regex_sets = _create_column_regex_all(self.regex_search_sets)
        not_search_regex = _create_column_regex_any(self.not_regex_search_strings)
        has_regex_search_strings = len(self.regex_search_strings) > 0
        has_regex_search_sets = len(self.regex_search_sets) > 0
        has_in_string_sets = len(self.in_string_sets) > 0
        has_no_not_regex_search_strings = len(self.not_regex_search_strings) == 0
        has_no_not_in_strings = len(self.not_in_strings) == 0
        columns_of_interest = [original_name for original_name, modified_name in name_map.items() if \
                                (
                                 (has_regex_search_strings and re.search(search_regex, modified_name)) or \
                                 (has_regex_search_sets and re.search(search_regex_sets, modified_name)) or \
                                 any([word in modified_name for word in self.in_strings]) or \
                                 (has_in_string_sets and any([all([word in modified_name for word in word_set]) for word_set in self.in_string_sets])) or \
                                 any([word == modified_name for word in self.exact_strings])
                                ) and \
                                (has_no_not_regex_search_strings or not re.search(not_search_regex, modified_name)) and \
                                (has_no_not_in_strings or all([word not in modified_name for word in self.not_in_strings]))]
        return columns_of_interest

class ValueMatcher():
    def __init__(self, values_type: None|str = None, values_regex: None|str = None, values_inverse_regex: None|str = None):
        self.values_type = values_type if isinstance(values_type, str) else ''
        self.values_regex = values_regex if isinstance(values_regex, str) else ''
        self.values_inverse_regex = values_inverse_regex if isinstance(values_inverse_regex, str) else ''
    
    def series_match(self, series: pandas.Series) -> pandas.Series:
        """Return a mask for the series based on type and regex matching.
        
        "values_regex" and "values_inverse_regex" are mutually exclusive and "values_regex" will take precedence if both are given. 
        "values_type" and one of the regex parameters can both be used, the intermediate selectors are ANDed together.
        
        Args:
            series: series to match values based on type and/or regex.
        
        Returns:
            A pandas Series the same length as "series" with Boolean values that can be used to select the matching values in the series.
        """
        if self.values_regex:
            regex_match = series.str.fullmatch(self.values_regex, na=False)
        elif self.values_inverse_regex:
            regex_match = ~series.str.fullmatch(self.values_inverse_regex, na=True)
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

class ColumnFinder:
    """Used to find columns in a DataFrame that match a NameMatcher and values in the column that match a ValueMatcher.
    """
    def __init__(self, standard_name: str, name_matcher: NameMatcher, value_matcher: ValueMatcher):
        self.standard_name = standard_name
        self.name_matcher = name_matcher
        self.value_matcher = value_matcher
    
    def name_dict_match(self, name_map):
        return self.name_matcher.dict_match(name_map)
    
    def values_series_match(self, series):
        return self.value_matcher.series_match(series)
    
    




def _find_matching_rows(series: pandas.Series, column_finder_list: list[ColumnFinder]) -> pandas.Series:
    """Return a mask for the series based on its list of matching attributes.
    
    The series should have been previously matched to some attributes for matching based on its 
    column name. It could be matched to more than one. If there are multiple matches, i.e. 
    there are multiple ColumnFinders in the column_finder_list, then the resulting masks 
    for each are ORed together. So a column name like "KEGG/PubChem" would match either 
    KEGG values or PubChem values.
    
    Args:
        series: series to match values.
        column_finder_list: a list of ColumnFinders.
        
    Returns:
        A pandas Series the same length as "series" with Boolean values that can be used to select the matching values in the series.
    """
    cumulative_mask = pandas.Series([False]*len(series), index=series.index)
    for column_finder in column_finder_list:
        mask = column_finder.values_series_match(series)
        cumulative_mask = cumulative_mask | mask
    
    return cumulative_mask


def _get_indexes_to_shift(column_to_shift_into: pandas.Series, rows_under_study: pandas.Series, 
                          column_finder_list: list[ColumnFinder], ignore_shift_values: bool = False) -> pandas.Index:
    """Get the indexes from rows_inder_study that can be shifted into the column_to_shift_into.
    
    rows_under_study have been identified as not appropriate for the column they are in, so test to see if they 
    are appropriate for the column_to_shift_into. By default, values will not be marked as shiftable if the 
    value currently in column_to_shift_into is appropriate for that column, but if ignore_shift_values is True, 
    then they will be ignored.
    
    Args:
        column_to_shift_into: the Series we are attempting to shift the values of rows_under_study into.
        rows_under_study: the current rows we are looking to shift into a different column in metabolites_df.
        column_finder_list: the list of ColumnFinders to be used to match values for the column_to_shift_into.
        ignore_shifted_values: if False, matching values in the column_to_shift_into won't allow values to be shifted into the column.
                               if True, ignore the values currently in the column_to_shift_into and shift onl based on values in rows_under_study.
    
    Returns:
        The indexes of rows_under_study that match the attributes in column_finder_list and therefore can be shifted.
    """
    rows_to_shift_into = column_to_shift_into.loc[rows_under_study.index]
    
    # See which values in the rows_under_study match against the column_to_shift_into's regexes.
    bad_values_match_in_shift_column = _find_matching_rows(rows_under_study, column_finder_list)
    # See if the values already in the indexes that want to shift match against the column_to_shift_into's regexes.
    shift_values_match_in_shift_column = _find_matching_rows(rows_to_shift_into, column_finder_list)
    
    if ignore_shift_values:
        rows_to_shift_mask = bad_values_match_in_shift_column
    else:
        rows_to_shift_mask = (bad_values_match_in_shift_column) & (~shift_values_match_in_shift_column)

    return rows_under_study[rows_to_shift_mask].index


def _update_temp_rows_to_shift(temp_rows_to_shift: dict, indexes: list, column_locations: list) -> dict:
    """Update temp_rows_to_shift with the indexes and column_locations.
    
    Just a simple utility to update the intermediate temp_rows_to_shift structure in this particular way.
    
    Args:
        temp_rows_to_shift: a dictionary of indexes to a list of column locations.
        indexes: the indexes to modify or add to the dict.
        column_locations: the column locations to assocaite with the indexes in the dict.
        
    Returns:
        The updated temp_rows_to_shift dict.
    """
    for index in indexes:
        if index in temp_rows_to_shift:
            temp_rows_to_shift[index] = list(set(temp_rows_to_shift[index] + column_locations))
        else:
            temp_rows_to_shift[index] = column_locations.copy()
    return temp_rows_to_shift


def _determine_rows_to_shift(metabolites_df: pandas.DataFrame, indexes_to_propagate: list[str|int], column_location: int, 
                             column_name_to_column_finders: dict[str, list[ColumnFinder]], 
                             cols_to_shift: dict, rows_to_shift: dict, 
                             right_shift: bool = True, last_unnamed: bool = False):
    """Given the column_location starting point and indexes_to_propagate, determine which rows can be shifted in metabolites_df.
    
    indexes_to_propagate have been determined that they need to be shifted in the column_location, so they 
    have to be tested to see if they can be. They can be shifted either because their neighbor also needs to 
    be shifted, or because their neighbor is a NA value. The testing must be propagated until a NA value is 
    found in the row, or a value is found that cannot be shifted, meaning the tested subsection of row cannot 
    be shifted. The last_unnamed parameter indicates that the column_location is the last one on the right and 
    has a name like "Unnamed: 0" or "". The assumption with columns like this are that all of their values need to be 
    shifted, so extra logic and checking is warranted to allow the shift.
    
    Args:
        metabolites_df: the DataFrame where rows are to be dhifted.
        indexes_to_propagate: the indexes in the DataFrame that have been determined need to be shifted from column_location.
        column_location: the location of the column in metabolites_df that has been identified as needing some values shifted.
        column_name_to_column_finders: a mapping of column names in metabolites_df to their list ColumnFinders.
        cols_to_shift: column locations are keys and the values are a list of indexes that have already been determined need to be shifted.
        ros_to_shift: indexes are keys and the values are a list of column locations that have already been determined need to be shifted.
        right_shift: if True, evaluate whether the indexes_to_propagate can be shifted right, else left.
        last_unnamed: if the column is in the last location and has a name like "Unnamed: 0" or "" then apply additional special case logic.
    
    Returns:
        A dicitonary where the keys are indexes and the values are a list of column locations where the indexes need to be shifted.
    """
    if right_shift:
        increment = 1
        loop_condition = lambda i: column_location + i + 1 < metabolites_df.shape[1]
        range_calculator = lambda i: list(range(column_location, column_location + i + 1))
    else:
        increment = -1
        loop_condition = lambda i: column_location + i > 1
        range_calculator = lambda i: list(range(column_location + i, column_location + 1))
    
    i = 0
    settled_indexes = []
    temp_rows_to_shift = {}
    while len(indexes_to_propagate) > 0 and loop_condition(i):
        next_column_name = metabolites_df.columns[column_location + i + increment]
        prev_column_name = metabolites_df.columns[column_location + i]
        bad_rows = metabolites_df.loc[indexes_to_propagate, prev_column_name]
        if column_name_to_column_finders[next_column_name]:
            indexes = _get_indexes_to_shift(metabolites_df.loc[:, next_column_name], bad_rows, 
                                            column_name_to_column_finders[next_column_name], 
                                            last_unnamed)
            NA_mask = metabolites_df.loc[:, next_column_name].isna()
            NA_indexes = metabolites_df.loc[:, next_column_name][NA_mask].index
            # If the column to the right has values that need to move left, 
            # and the column to the left does NOT have values that need to move right, 
            # then right_indexes can be settled. i must be 0, otherwise we know
            # already that there are values in the left column wanting right. 
            # This also applies when shifting left.
            if i == 0 and not last_unnamed:
                series = metabolites_df.loc[:, next_column_name]
                column_mask = _find_matching_rows(series, column_name_to_column_finders[next_column_name]) | series.isna()
                next_bad_rows = series[~column_mask]
                indexes_needing_opposite_shift = _get_indexes_to_shift(metabolites_df.loc[:, prev_column_name], 
                                                                       next_bad_rows, 
                                                                       column_name_to_column_finders[prev_column_name])
                indexes_needing_shift = []
                if column_location + i - increment in cols_to_shift:
                    indexes_needing_shift = cols_to_shift[column_location + i - increment]
                settled_indexes = [index for index in indexes if (index not in indexes_needing_shift and \
                                                                 index in indexes_needing_opposite_shift) or \
                                                                 index in NA_indexes]
            else:
                settled_indexes = [index for index in indexes if index in NA_indexes or \
                                                                          (index in rows_to_shift and \
                                                                          column_location + i - 1 in rows_to_shift[index])]
            indexes_to_propagate = [index for index in indexes if index not in settled_indexes]
            temp_rows_to_shift = _update_temp_rows_to_shift(temp_rows_to_shift, settled_indexes, range_calculator(i))
            
        i += increment
        
    # Add a special case for AN004201. If the propagation goes all the way to the Metabolite column and 
    # the value in that column is the only number in the column, overwrite the value.
    if last_unnamed and column_location + i == 1 and len(indexes_to_propagate) == 1:
        series = metabolites_df.iloc[:, 0]
        mask = ColumnFinder('', NameMatcher(), ValueMatcher(values_type = 'non-numeric')).values_series_match(series)
        numeric_metabolites = series[~mask]
        if numeric_metabolites.shape[0] == 1 and \
           (str(indexes_to_propagate[0]) in numeric_metabolites.index or \
            indexes_to_propagate[0] in numeric_metabolites.index):
            temp_rows_to_shift = _update_temp_rows_to_shift(temp_rows_to_shift, indexes_to_propagate, list(range(0, column_location + 1)))
    
    return temp_rows_to_shift





def _update_rows_and_cols_to_shift(rows_to_shift: dict, cols_to_shift: dict, temp_rows_to_shift: dict) -> Tuple[dict,dict]:
    """Update rows_to_shift and cols_to_shift with the values in temp_rows_to_shift.
    
    Just a simple utility to update the persistent structures from the intermediate 
    temp_rows_to_shift structure in this particular way.
    
    Args:
        rows_to_shift: a dictionary of indexes to a list of column locations.
        cols_to_shift: column locations are keys and the values are a list of indexes. 
                       The same information as in rows_to_shift, but from a column persepective.
        temp_rows_to_shift: a dictionary of indexes to a list of column locations.
        
    Returns:
        The updated rows_to_shift and cols_to_shift dicts.
    """
    for index, column_locs in temp_rows_to_shift.items():
        if index in rows_to_shift:
            rows_to_shift[index] += column_locs
        else:
            rows_to_shift[index] = column_locs
        for column_loc in column_locs:
            if column_loc in cols_to_shift:
                cols_to_shift[column_loc].add(index)
            else:
                cols_to_shift[column_loc] = {index}
    return rows_to_shift, cols_to_shift


def _determine_rows_and_cols_to_shift(metabolites_df: pandas.DataFrame, column_finders: dict[str, ColumnFinder]) -> Tuple[dict, dict, dict, dict]:
    """Go through each column of metabolites_df and look for values to shift within a row.
    
    The main work horse of fixing row shifted values in a METABOLITES data block. Uses 
    column_finders to match column names and values in the columns to be shifted within a row.
    
    Args:
        metabolites_df: the METABOLITES data block of a mwTab file in DataFrame form.
        column_finders: the keys are generalized column names and the values are ColumnFinders.
        
    Returns:
        rows_to_shift_right, rows_to_shift_left, cols_to_shift_right, cols_to_shift_left
        rows_to_shift are dicts with keys for each index to shift and a list of column locations as values.
        cols_to_shift are dicts with keys for each column location and a list of indexes as values.
        Both are the same information just in 2 different perspectives for easier searching downstream.
    """
    column_name_to_column_finders = {column:[] for column in metabolites_df.columns}
    column_name_map = {column:column.lower() for column in metabolites_df.columns}
    for column_finder in column_finders.values():
        columns_of_interest = column_finder.name_dict_match(column_name_map)
        for column in columns_of_interest:
            column_name_to_column_finders[column].append(column_finder)
    
    if any(len(column_name_to_column_finders[column]) > 0 for column in metabolites_df.columns[1:]):
        rows_to_shift_right = {}
        rows_to_shift_left = {}
        cols_to_shift_right = {}
        cols_to_shift_left = {}
        for column_name in metabolites_df.columns[1:]:
            column_location = metabolites_df.columns.get_loc(column_name)
            series = metabolites_df.loc[:, column_name]
            temp_rows_to_shift_right = {}
            temp_rows_to_shift_left = {}
            
            if column_name_to_column_finders[column_name]:
                matching_row_mask = _find_matching_rows(series, column_name_to_column_finders[column_name]) | series.isna()
                if not matching_row_mask.all():
                    bad_rows = series[~matching_row_mask]
                    
                    if column_location < metabolites_df.shape[1] - 1:
                        temp_rows_to_shift_right = _determine_rows_to_shift(metabolites_df, list(bad_rows.index), column_location, 
                                                                            column_name_to_column_finders, 
                                                                            cols_to_shift_right, rows_to_shift_right, 
                                                                            right_shift = True, last_unnamed = False)
                    
                    if column_location > 1:
                        temp_rows_to_shift_left = _determine_rows_to_shift(metabolites_df, list(bad_rows.index), column_location, 
                                                                           column_name_to_column_finders, 
                                                                           cols_to_shift_left, rows_to_shift_left, 
                                                                           right_shift = False, last_unnamed = False)
            elif ("Unnamed" in column_name or column_name == "" or re.match(r'\{\{\{_\d+_\}\}\}', column)) \
                 and column_location == metabolites_df.shape[1] - 1:
                trimmed_series = series[~series.isna()]
                # Some unnamed are legit columns just missing the name, detect this by requiring most of the values to be NA.
                if len(trimmed_series)/len(series) > .5:
                    continue
                
                temp_rows_to_shift_left = _determine_rows_to_shift(metabolites_df, list(trimmed_series.index), column_location, 
                                                                   column_name_to_column_finders, 
                                                                   cols_to_shift_left, rows_to_shift_left, 
                                                                   right_shift = False, last_unnamed = True)
            else:
                continue
                        
            # To try and eliminate an issue where the same value could be shifted left or right, 
            # we compare the non-bad values in left and right and give favor to the side that 
            # has a more restrictive regex.
            a_right_column_exists = column_location + 1 < metabolites_df.shape[1] and \
                                    column_name_to_column_finders[metabolites_df.columns[column_location + 1]]
            a_left_column_exists = column_location - 1 > 1 and \
                                   column_name_to_column_finders[metabolites_df.columns[column_location - 1]]
            if a_right_column_exists and a_left_column_exists:
                right_column_name = metabolites_df.columns[column_location + 1]
                left_column_name = metabolites_df.columns[column_location - 1]
                
                right_column_not_bad_rows = metabolites_df.loc[:, right_column_name][matching_row_mask]
                left_column_not_bad_rows = metabolites_df.loc[:, left_column_name][matching_row_mask]
                
                right_values_left_match = _find_matching_rows(right_column_not_bad_rows, column_name_to_column_finders[left_column_name]).any()
                left_values_right_match = _find_matching_rows(left_column_not_bad_rows, column_name_to_column_finders[right_column_name]).any()
                if right_values_left_match and not left_values_right_match:
                    temp_rows_to_shift_left = {index:value for index, value in temp_rows_to_shift_left.items() if index not in temp_rows_to_shift_right}

                if left_values_right_match and not right_values_left_match:
                    temp_rows_to_shift_right = {index:value for index, value in temp_rows_to_shift_right.items() if index not in temp_rows_to_shift_left}
            
            
            rows_to_shift_right, cols_to_shift_right = _update_rows_and_cols_to_shift(rows_to_shift_right, 
                                                                                      cols_to_shift_right, 
                                                                                      temp_rows_to_shift_right)
            rows_to_shift_left, cols_to_shift_left = _update_rows_and_cols_to_shift(rows_to_shift_left, 
                                                                                    cols_to_shift_left, 
                                                                                    temp_rows_to_shift_left)
            
        # Keep only unique values and sort them.
        rows_to_shift_right = {int(index):sorted(list(set(value))) for index, value in rows_to_shift_right.items()}
        rows_to_shift_left = {int(index):sorted(list(set(value))) for index, value in rows_to_shift_left.items()}
        
        return rows_to_shift_right, rows_to_shift_left, cols_to_shift_right, cols_to_shift_left
    return {}, {}, {}, {}







def _find_dominoe_columns(columns_of_interest: list[int], all_columns: list[int], leftwise: bool = True) -> list[int]:
    """Find all neighboring columns to columns in columns_of_interest from all_columns.
    
    columns_of_interest can't be shifted, so find all other columns in all_columns that 
    also can't be shifted due to their proximity to the columns_of_interest. Note that 
    there can be repeated values in the returned list, and they are in no certain order.
    
    Args:
        columns_of_interest: the columns to look for neighbors to.
        all_columns: the full set of shiftable columns.
        leftwise: if True, assume a left shit, otherwise assume a right shift.
    
    Returns:
        The list of columns that cannot be shifted, including columns_of_interest.
    """
    increment = 1
    if not leftwise:
        increment = -1
        
    dominoe_columns = columns_of_interest.copy()
    for column in columns_of_interest:
        next_column = column + increment
        while next_column in all_columns:
            dominoe_columns.append(next_column)
            next_column = next_column + increment
    return dominoe_columns


def _modify_shifted_df(metabolites_df: pandas.DataFrame, shifted_metabolites_df: pandas.DataFrame, 
                       rows_to_shift_right: dict, rows_to_shift_left: dict, 
                       cols_to_shift_right: dict, cols_to_shift_left: dict, 
                       shift_right: bool = True) -> Tuple[pandas.DataFrame, list[str]]:
    """Modify metabolites_df to do the necessary shifts according to the given shift parameters.
    
    metabolites_df should be the unchanged METABOLITES block DataFrame, and shifted_metabolites_df 
    should be a shifted version of it. The shift should match the state of shift_right. So if 
    shift_right is True, then shifted_metabolites_df should be a version where all rows are shifted 
    1 space to the right. This function does the shifting, but also settles some last collisions 
    that can happen between values needing to shift left and right.
    
    Args:
        metabolites_df: the unchanged METABOLITES block DataFrame.
        shifted_metabolites_df: the original METABOLITES block DataFrame shifted left or right.
        rows_to_shift_right: a dictionary of indexes to a list of column locations for rows needing to be shifted to the right.
        rows_to_shift_left: a dictionary of indexes to a list of column locations for rows needing to be shifted to the left.
        cols_to_shift_right: column locations are keys and the values are a list of indexes.
                             The same information as in rows_to_shift_right, but from a column persepective.
        cols_to_shift_left: column locations are keys and the values are a list of indexes.
                            The same information as in rows_to_shift_left, but from a column persepective.
        shift_right: if True, try to shift the indexes in rows_to_shift_right to the right, 
                     otherwise try to shift the indexes in rows_to_shift_left to the left
    
    Returns:
        metabolites_df modified and a list of all messages about what was shifted.
    """
    
    if shift_right:
        rows_to_shift_primary = rows_to_shift_right
        rows_to_shift_secondary = rows_to_shift_left
        cols_to_shift_primary = cols_to_shift_right
        cols_to_shift_secondary = cols_to_shift_left
        message_string = 'Right'
        dominoe_direction = False
        increment = 1
    else:
        rows_to_shift_primary = rows_to_shift_left
        rows_to_shift_secondary = rows_to_shift_right
        cols_to_shift_primary = cols_to_shift_left
        cols_to_shift_secondary = cols_to_shift_right
        message_string = 'Left'
        dominoe_direction = True
        increment = -1

    messages = []
    for index, column_locs in rows_to_shift_primary.items():
        if index in rows_to_shift_secondary:
            primary_column_coverage = [column_loc + 1 for column_loc in column_locs]
            secondary_column_coverage = [column_loc - 1 for column_loc in rows_to_shift_secondary[index]]
            secondary_competing_columns = [column_loc for column_loc in secondary_column_coverage if column_loc in primary_column_coverage]
            secondary_competing_columns += [column_loc for column_loc in rows_to_shift_secondary[index] if column_loc in column_locs]
            secondary_competing_columns = sorted(secondary_competing_columns)
            # If there are competing shifts, give preference to the side with the most columns.
            if secondary_competing_columns:
                if len(column_locs) > len(rows_to_shift_secondary[index]):
                    # If a column can't be left/right shifted then the ones immediately to the right/left also can't be left/right shifted.
                    secondary_dominoe_columns = _find_dominoe_columns(secondary_competing_columns, rows_to_shift_secondary[index], not dominoe_direction)
                    rows_to_shift_secondary[index] = [column_loc for column_loc in rows_to_shift_secondary[index] if column_loc not in secondary_dominoe_columns]
                elif len(rows_to_shift_secondary[index]) > len(column_locs):
                    primary_competing_columns = [column_loc for column_loc in primary_column_coverage if column_loc in secondary_column_coverage]
                    primary_dominoe_columns = _find_dominoe_columns(primary_competing_columns, rows_to_shift_primary[index], dominoe_direction)
                    rows_to_shift_primary[index] = [column_loc for column_loc in rows_to_shift_primary[index] if column_loc not in primary_dominoe_columns]
                else:
                    sum_of_primary_rows_in_column = sum([len(cols_to_shift_primary[column]) for column in rows_to_shift_primary[index]])
                    sum_of_secondary_rows_in_column = sum([len(cols_to_shift_secondary[column]) for column in rows_to_shift_secondary[index]])
                    if sum_of_primary_rows_in_column > sum_of_secondary_rows_in_column:
                        secondary_dominoe_columns = _find_dominoe_columns(secondary_competing_columns, rows_to_shift_secondary[index], not dominoe_direction)
                        rows_to_shift_secondary[index] = [column_loc for column_loc in rows_to_shift_secondary[index] if column_loc not in secondary_dominoe_columns]
                    elif sum_of_secondary_rows_in_column > sum_of_primary_rows_in_column:
                        primary_competing_columns = [column_loc for column_loc in primary_column_coverage if column_loc in secondary_column_coverage]
                        primary_dominoe_columns = _find_dominoe_columns(primary_competing_columns, rows_to_shift_primary[index], dominoe_direction)
                        rows_to_shift_primary[index] = [column_loc for column_loc in rows_to_shift_primary[index] if column_loc not in primary_dominoe_columns]
                    
        filtered_column_locs = rows_to_shift_primary[index]
        for i, column_loc in enumerate(filtered_column_locs):
            metabolites_df.iloc[index, column_loc + increment] = shifted_metabolites_df.iloc[index, column_loc + increment]
            if shift_right: 
                if i == 0 or column_loc - 1 != filtered_column_locs[i-1]:
                    metabolites_df.iloc[index, column_loc] = pandas.NA
            else:
                if not (index in rows_to_shift_secondary and column_loc - 1 in rows_to_shift_secondary[index]) and \
                   (i+1 == len(filtered_column_locs) or column_loc - 1 != filtered_column_locs[i-1]):
                       metabolites_df.iloc[index, column_loc] = pandas.NA
        if filtered_column_locs:
            messages.append(f"{message_string} shift done on row {index}.")
    
    return metabolites_df, messages

def fix_row_shift(metabolites_df: pandas.DataFrame, column_finders: dict[str, ColumnFinder]) -> Tuple[pandas.DataFrame, list[str]]:
    """Find and fix and rows in metabolites_df that need shifting left or right.
    
    metabolites_df should be the METABOLITES block of an mwTab file in DataFrame form. 
    column_finders is a dictionary of column names with values that are ColumnFinders.
    
    Args:
        metablites_df: the METABOLITES block of an mwTab file in DataFrame form.
        column_finders: a dict where the keys are column names and the values are ColumnFinders.
        
    Returns:
        metabolites_df modified if needed and a list of all messages about what was shifted. 
        If the messages are an empty list, then there should be no change to metabolites_df.
    """
    rows_and_columns_to_shift = _determine_rows_and_cols_to_shift(metabolites_df, column_finders)
    rows_to_shift_right, rows_to_shift_left, cols_to_shift_right, cols_to_shift_left = rows_and_columns_to_shift
    
    messages = []
    if rows_to_shift_right or rows_to_shift_left:
        left_shift_metabolites_df = metabolites_df.shift(-1, axis=1).copy()
        right_shift_metabolites_df = metabolites_df.shift(1, axis=1).copy()
        
        metabolites_df, modify_messages = _modify_shifted_df(metabolites_df, right_shift_metabolites_df, 
                                                             rows_to_shift_right, rows_to_shift_left, 
                                                             cols_to_shift_right, cols_to_shift_left, 
                                                             shift_right = True)
        messages += modify_messages
        metabolites_df, modify_messages = _modify_shifted_df(metabolites_df, left_shift_metabolites_df, 
                                                             rows_to_shift_right, rows_to_shift_left, 
                                                             cols_to_shift_right, cols_to_shift_left, 
                                                             shift_right = False)
        messages += modify_messages
        
        # After shifting, some columns on the end will now be all NA with no name, so remove them.
        no_name_columns = [column for column in metabolites_df.columns if column == '' or re.match(r'\{\{\{_\d+_\}\}\}', column)]
        if no_name_columns:
            no_name_df = metabolites_df.loc[:, no_name_columns]
            null_columns = no_name_df.isna().all()
            columns_to_drop = null_columns[null_columns].index
            if len(columns_to_drop) > 0:
                metabolites_df = metabolites_df.drop(columns = columns_to_drop)
        
        
    return metabolites_df, messages


