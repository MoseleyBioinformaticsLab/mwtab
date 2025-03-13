# -*- coding: utf-8 -*-
"""
The repair command.
"""

import re
import traceback
import sys
from itertools import combinations
import html
import warnings
warnings.filterwarnings("ignore", module="fuzzywuzzy")

import pandas
import numpy
import fuzzywuzzy.fuzz


from . import mwtab
from .validator import METABOLITES_REGEXS
from .duplicates_dict import DuplicatesDict



VERBOSE = False
LIST_SPLIT_REGEX = r'/\s*|\s*and\s*|;\s*|,\s*|//\s*|\s*\&\s*'
# TODO replace these NA values with the empty string.
# Sometimes 0 is also an NA value, but it can be hard tell unless you see it in an ID column like KEGG ID or something.
# Consider looking for 0 in non-numeric rows and then replacing with the empty string.
# Note the slightly different hyphen character. It is not a duplicate.
NA_VALUES = ['', '-', '−', '--', '---', '.', ',',
             'NA', 'na', 'n.a.', 'N.A.', 'n/a', 'N/A', '#N/A', 'NaN', 'nan', 'N', 'null', 'Null', 'NULL', 'NF',
             'No result', 'NOT Found in Database', 'No ID', 'no data', 'unknown', 'undefined', 'No record', 'NIDB',
             'Not available', 'TBC', 'Internal Standard', 'Intstd', 'internal standard', 'Internal standard',
             'Spiked Stable Isotope Labeled Internal Standards', 'Int Std']

# AN000198 has the bottom 4 rows shifted after the ri_type column, need to be shifted right.
# TODO when trying to look for columns to change to a common name, probably filter out ones with names in parantheses, 'RSD (Mass, ppm)' AN002782
# TODO AN004201 has a row in METABOLITES that is shifted to the right and has a nonsense value for the metabolite, fix this.
# AN004512 has some values that are shifted to the right as well.
# AN003295 has the kegg and pubchem columns swapped, fix this.
# When renaming metabolites columns to standard names, make sure there is only 1 match. Think about what to do for multiple matches.
# Create a documentation page that gives examples of all of the known/standardized metabolites columns.
# Can use literal_eval to read in lists from strings.
# from ast import literal_eval      temp = literal_eval("['a', 'b']")
# Add repair block to mwtab file that details changes made.
# AN000370 Part of the Metabolites column is duplicated in the column next to it. 
#    Look for columns that are duplicates or partial duplicates of Metabolites and set them to NA. 
#    Just do it and look at the ones that get modified to see exactly what happens.
# AN000084 has another table pasted to the bottom in METABOLITES and it's just wrong. Not sure how to fix or if it should be.
# AN003295 has a KEGG and pubchem column, but they are mislabled. The KEGG column has pubchem values and vice versa.

    
WRAP_STRING = r'[^a-zA-Z0-9]'
def create_column_regex_any(match_strings):
    """
    """
    regex = '|'.join(['(' + WRAP_STRING + '|^)' + match_string + '(' + WRAP_STRING + '|$)' for match_string in match_strings])
    return regex

def create_column_regex_all(match_strings_list):
    """
    """
    regex = '|'.join([''.join(['(?=.*(' + WRAP_STRING + '|^)' + match_string + '(' + WRAP_STRING + '|$))' for match_string in match_strings]) for match_strings in match_strings_list])
    return regex

def column_name_matching(column_names, lowered_column_names, 
                         regex_search_strings, not_regex_search_strings, regex_search_sets,
                         in_strings, not_in_strings, in_string_sets, exact_strings):
    """
    """
    search_regex = create_column_regex_any(regex_search_strings)
    search_regex_sets = create_column_regex_all(regex_search_sets)
    not_search_regex = create_column_regex_any(not_regex_search_strings)
    has_regex_search_strings = len(regex_search_strings) > 0
    has_regex_search_sets = len(regex_search_sets) > 0
    has_in_string_sets = len(in_string_sets) > 0
    has_no_not_regex_search_strings = len(not_regex_search_strings) == 0
    has_no_not_in_strings = len(not_in_strings) == 0
    columns_of_interest = [column_names[i] for i, column_name in enumerate(lowered_column_names) if \
                            (
                             (has_regex_search_strings and re.search(search_regex, column_name)) or \
                             (has_regex_search_sets and re.search(search_regex_sets, column_name)) or \
                             any([word in column_name for word in in_strings]) or \
                             (has_in_string_sets and any([all([word in column_name for word in word_set]) for word_set in in_string_sets])) or \
                             any([word == column_name for word in exact_strings])
                            ) and \
                            (has_no_not_regex_search_strings or not re.search(not_search_regex, column_name)) and \
                            (has_no_not_in_strings or all([word not in column_name for word in not_in_strings]))]
    return columns_of_interest


{
 "Duplicate keys in SSF Additional sample data, and order matters." : [],
 "Internal representations do not match." : ["partially repairable"], # This is more generalized and some can be repaired.
 "WORKBENCH header has duplicate keys." : ["partially repairable"], # If the key also has duplicate values can be repaired. For different values just add a prefix maybe?
 "Non key value element(s) in WORKBENCH header missing." : [], # This is more of an issue of having non key value elements at all and whether that's an error or not.
 # TODO change parsing so it pairs up the filename before the datatrack_ID and give the filename the key of "DATATRACK_ID:###". 
 # In other words, make the key the associated datatrack_ID key value pair.
 # Need to think about this and look at it some more. Would be much easier to just ignore this stuff. Not sure how much value it has.
 # TODO think about adding a more complicated check for data where you compare the type of the data in each column, and size, etc
 # and use that information to make corrections. This could work for AN000593 (search for this in temp.py it explains the situation).
 # I added a type check that covers AN000593, but not sure if that's enough.
 "Duplicate sample names." : [],
 # AN000144 has the samples in SSF in a different order than what is in the header. The duplicated sample is near the middle of the list 
 # in the header, S00011727, but the missing sample is S00011706 which should be near the beginning, so it could be that all the names 
 # should be shifted so that S00011706 is in the correct spot, or it could be that S00011706 should replace one of the S00011727's. 
 # I don't think we can say, so I don't think this can be repaired.
 
 # TODO Duplicate factor names in SSF. AN000379   This looks handled. The duplicate names are prefixed after repair.
 # TODO Changed the tokenizer so factors split on ' | ' instead of '|'. Make sure that files still parse okay. This was to fix an issue in AN001296 where a factor had a pipe in it.
 # I have checked this on some files and it's fine, but haven't rigorously looked through them all or anything.
 }


def printing(message):
    """
    """
    if VERBOSE:
        print(message)


def data_block_element_fix(start_index, end_index, line, index, lines):
    """
    """
    if start_index and not end_index and re.match(r'.*"[^"]*\t[^"]*".*', line):
        tokens = []
        line_copy = line
        while '\t' in line_copy:
            tab_index = line_copy.index('\t')
            token = line_copy[0:tab_index]
            if token.count('"') % 2 == 0:
                tokens.append(token)
                line_copy = line_copy[tab_index+1:]
            else:
                line_copy = line_copy[0:tab_index] + line_copy[tab_index+1:]
        tokens.append(line_copy)
        
        lines[index] = '\t'.join([token.strip('" ') for token in tokens])
        return True
    return False


def is_index_in_range(index, low_end, high_end):
    if not low_end or not high_end:
        return False
    if index > low_end and index < high_end:
        return True
    return False

def is_in_a_metabolite_block(index, data_low, data_high, extended_low, extended_high, metab_low, metab_high, binned_low, binned_high):
    if is_index_in_range(index, data_low, data_high) or\
       is_index_in_range(index, extended_low, extended_high) or\
       is_index_in_range(index, metab_low, metab_high) or\
       is_index_in_range(index, binned_low, binned_high):
        return True

def line_should_have_two_letter_code(line, index, data_low, data_high, extended_low, extended_high, metab_low, metab_high, binned_low, binned_high):
    return (not line.startswith("#") and \
            not line.startswith("VERSION") and \
            not line.startswith("CREATED_ON") and \
            not line.startswith("SUBJECT_SAMPLE_FACTORS") and \
            not "METABOLITE_DATA" in line and \
            not "BINNED_DATA" in line and \
            not line == "" and \
            not is_in_a_metabolite_block(index, 
                                         data_low, data_high, 
                                         extended_low, extended_high, 
                                         metab_low, metab_high,
                                         binned_low, binned_high))

def line_is_numbers(line):
    tokens = line.split('\t')[1:]
    bool_list = []
    for token in tokens:
        try:
            float(token)
            bool_list.append(True)
        except Exception:
            if token == '' or token == 'NA':
                bool_list.append(True)
            else:
                bool_list.append(False)
    return all(bool_list)

def is_column_numeric(column):
    """
    """
    stripped_column = column.str.strip()
    stripped_column = stripped_column.str.strip('\u200e')
    old_NAs = (stripped_column.isna()) | (stripped_column.isin(NA_VALUES))
    column_to_numeric = pandas.to_numeric(stripped_column, errors='coerce')
    column_to_numeric_NAs = column_to_numeric.isna()
    column_to_numeric_nums = ~column_to_numeric_NAs
    new_NAs = column_to_numeric_NAs ^ old_NAs
    return column_to_numeric_nums.any() and len(stripped_column[new_NAs]) == 0

def is_column_ID(column):
    """
    Want to find IDs. For example, C17456.
    Patterns are, empty string, 0, a single letter followed by at least 4 numbers, and at least 4 numbers.
    """
    stripped_column = column.str.strip()
    stripped_column = stripped_column.str.strip('\u200e')
    return stripped_column.str.fullmatch(r'([a-zA-Z]\d{4,})|(\d{4,})|0|').all()
 
def column_has_kegg_ID(column):
    """
    Want to find KEGG IDs. For example, C17456.
    Pattern is a 'c' or 'C' followed by 5 numbers.
    """
    stripped_column = column.str.strip()
    stripped_column = stripped_column.str.strip('\u200e')
    return stripped_column.str.fullmatch(r'[cC]\d{5}').any()

def column_has_inchi_key(column):
    """
    Want to find InChIKeys. For example, UWPWWENWLZPQGU-WRFRXMDISA-O.
    """
    stripped_column = column.str.strip()
    stripped_column = stripped_column.str.strip('\u200e')
    return stripped_column.str.fullmatch(r'[a-zA-Z]{14}-[a-zA-Z]{10}-[a-zA-Z]').any()

def consolidate_values_to_column(data_df, table_name, match_words, element_regex):
    """
    """
    stripped_columns = {column:data_df.loc[:, column].str.strip() for column in data_df.columns}
    stripped_columns = {column:data.str.strip('\u200e') for column, data in stripped_columns.items()}
    # Look for elementes that aren't in the matching column and move them.
    matching_columns = [column for column in data_df.columns if all([word in column.lower() for word in match_words])]
    if matching_columns:
        matching_column = matching_columns[0]
        boolean_columns = {column:data.str.fullmatch(element_regex) for column, data in stripped_columns.items() if data.str.fullmatch(element_regex).any()}
        if len(boolean_columns) > 1:
            boolean_df = pandas.DataFrame(boolean_columns)
            rows_with_single_match = (boolean_df.sum(axis=1) == 1)
            # Add elements to the matching_column and set them to the empty string in the other columns.
            for column, data in boolean_columns.items():
                if column != matching_column:
                    relevant_rows = rows_with_single_match & data
                    data_df.loc[relevant_rows, matching_column] = data_df.loc[relevant_rows, column]
                    data_df.loc[relevant_rows, column] = ''
            print("Elements moved in '" + table_name + " to column: " + matching_column)

def string_to_list(value):
    """
    """
    stripped_value = value.strip()
    if stripped_value[0] == '[' and stripped_value[-1] == ']':
        stripped_value = stripped_value[1:-1]
    split_values = re.split(LIST_SPLIT_REGEX, stripped_value)
    return [value.strip("'\"") for value in split_values]

def are_list_values_nums(values):
    """
    """
    values = [value for value in values if not value in NA_VALUES and not pandas.isna(value)]
    return not any(pandas.isna(pandas.to_numeric(values, errors='coerce')))

def fix_missing_tab_in_table(data_df, table_name):
    """
    """    
    for i, column in enumerate([column for column in data_df.columns if column != 'Metabolite']):
        temp_column = data_df.loc[:, column].astype(str)
        value_counts = temp_column.value_counts(dropna=False)
        
        # If there is exactly 1 value that is different it is likely that 
        # the tabs in the row were off, so try to determine that and fix.
        condition = value_counts == 1
        if any(condition) and not all(condition) and len(value_counts) == 2:
            one_value = value_counts[condition].index[0]
            fix_value = value_counts[~condition].index[0]
            one_row_index = temp_column[temp_column == one_value].index[0]
            # If the one value is in the value that is there we will assume a tab is missing.
            if one_value != "" and fix_value != "" and one_value.startswith(fix_value):
                print("Missing tab fixed in '" + table_name + "'. Column: " + column)
                # Have to get the fix_value in it's natural data type.
                value_counts = data_df.loc[:, column].value_counts(dropna=False)
                new_fix_value = value_counts[value_counts != 1].index[0]
                data_df.loc[one_row_index, column] = new_fix_value
                # Need to pull out the value to fix the column to the right.
                second_fix_value = one_value.replace(fix_value, '', 1)
                # Fix all columns to the right.
                values_to_replace = list(data_df.loc[one_row_index, :].values)
                # We skipped the Metabolite column, so have to make the counts line up right.
                if i + 2 < data_df.shape[1]:
                    right_column = data_df.columns[i + 2]
                    # Try to find the type of the column to the right by looking above and below the index.
                    if one_row_index - 1 > 0:
                        good_index = one_row_index - 1
                    else:
                        good_index = one_row_index + 1
                    second_fix_value = type(data_df.loc[good_index, right_column])(second_fix_value)
                    data_df.loc[one_row_index, right_column] = second_fix_value
                    j = i + 3
                    while j < data_df.shape[1]:
                        right_column = data_df.columns[j]
                        data_df.loc[one_row_index, right_column] = type(data_df.loc[good_index, right_column])(values_to_replace[j-1])
                        j += 1
        
        # Look for bad spacing by looking for different types.
        temp_column = data_df.loc[:, column].astype(str)
        # Moved stripping into clean_columns
        # stripped_column = temp_column.str.strip()
        # stripped_column = stripped_column.str.strip('\u200e')
        old_NAs = (temp_column.isna()) | (temp_column.isin(NA_VALUES))
        column_to_numeric = pandas.to_numeric(temp_column, errors='coerce')
        column_to_numeric_NAs = column_to_numeric.isna()
        column_to_numeric_nums = ~column_to_numeric_NAs
        new_NAs = column_to_numeric_NAs ^ old_NAs
        NAs_removed = temp_column[~new_NAs]
        number_of_new_NAs = len(temp_column[new_NAs])
        # Must have at least 1 numeric value and there can only be 1 non-numeric value.
        # if column_to_numeric_nums.any() and len(temp_column[new_NAs]) == 1:
        if column_to_numeric_nums.sum() > number_of_new_NAs and number_of_new_NAs > 0:
            for bad_index, bad_value in temp_column[new_NAs].items():
                # If the bad_value is in the metabolite then ignore it.
                if bad_value in data_df.loc[bad_index, 'Metabolite']:
                    continue
                # Find the longest substring that can be converted to numeric to aide in splitting the value.
                # If the value is a split value, ex 16857/15700 and has a corresponding metabolite name, threonine/homoserine, then ignore.
                # Also ignore if the column name ends in an 's'. The idea is that a plural could indicate multiple values.
                list_check = string_to_list(bad_value)
                # TODO maybe see if any value in the row can be split as well, not just metabolite. Example, AN001792 has a 'Species' column that lines up with m/z values.
                metabolite_check = re.split(LIST_SPLIT_REGEX, data_df.loc[bad_index, 'Metabolite'])
                list_values_are_nums = are_list_values_nums(list_check)
                if (len(list_check) == len(metabolite_check) and list_values_are_nums) or\
                   (list_values_are_nums and len(list_check) > 1 and column.lower().strip().endswith('s')):
                    continue
                # If the bad value is a range, i.e. separated by a hyphen or underscore then skip.
                range_check = re.split(r'\s*-\s*|\s*_\s*', bad_value)
                if len(range_check) == 2 and not any(pandas.isna(pandas.to_numeric(range_check, errors='coerce'))):
                    continue
                length = len(bad_value) + 1
                all_substrings = pandas.Series([bad_value[x:y] for x, y in combinations(range(length), r=2)])
                subs_to_numeric = pandas.to_numeric(all_substrings, errors='coerce')
                convertable_subs = all_substrings[~subs_to_numeric.isna()]
                if len(convertable_subs) == 0:
                    continue
                good_value = convertable_subs.loc[convertable_subs.str.len().idxmax()]
                # If the good_value does not have the same numeric type as the column then remove the whole value instead.
                if NAs_removed.str.contains('.', regex=False).all() and not '.' in good_value or \
                   not NAs_removed.str.contains('.', regex=False).all() and '.' in good_value:
                    data_df.iloc[bad_index, i+1] = numpy.nan
                    print("Non-numeric value removed from column " + column + " in " + table_name)
                    continue
                # Determine which side has the extra characters.
                re_match = re.match(r'(.*)' + good_value + r'(.*)', bad_value)
                left_match = re_match.group(1)
                right_match = re_match.group(2)
                if left_match and right_match:
                    continue
                elif left_match:
                    if i > 0 and not is_column_numeric(data_df.iloc[:, i]) and not is_column_ID(data_df.iloc[:, i]):
                        data_df.iloc[bad_index, i] = data_df.iloc[bad_index, i] + left_match
                        print("String moved to column '" + data_df.columns[i] + " in " + table_name)
                    data_df.iloc[bad_index, i+1] = good_value.strip()
                    print("Numeric type separated in '" + table_name + "'. Column: " + column)
                elif right_match:
                    # A question mark at the end of the value can indicate that the value is uncertain. Limit this to ID columns.
                    if right_match.strip() == '?' and is_column_ID(temp_column[column_to_numeric_nums]):
                        continue
                    if i + 2 < data_df.shape[1] and not is_column_numeric(data_df.iloc[:, i+2]) and not is_column_ID(data_df.iloc[:, i+2]):
                        data_df.iloc[bad_index, i+2] = right_match + data_df.iloc[bad_index, i+2]
                        print("String moved to column '" + data_df.columns[i+2] + " in " + table_name)
                    data_df.iloc[bad_index, i+1] = good_value.strip()
                    print("Numeric type separated in '" + table_name + "'. Column: " + column)
    
    return data_df

def compute_fuzz_ratios(list1, list2, inclusion_ratio=90):
    """
    """
    fuzz_ratios = {}
    for element1 in list1:
        temp_dict = {}
        lowered_element1 = element1.lower()
        for element2 in list2:
            
            lowered_element2 = element2.lower()
            if len(element1) == len(element2):
                replacement1 = lowered_element1.replace('^', "'")
                replacement1 = replacement1.replace('_', ',')
                replacement2 = lowered_element1.replace("'", "^")
                replacement2 = replacement2.replace(',', '_')
                replacement3 = lowered_element1.replace("'", "^")
                replacement3 = replacement3.replace('_', ',')
                replacement4 = lowered_element1.replace("^", "'")
                replacement4 = replacement4.replace(',', '_')
                if replacement1 == lowered_element2 or \
                   replacement2 == lowered_element2 or \
                   replacement3 == lowered_element2 or \
                   replacement4 == lowered_element2:
                    temp_dict[element2] = 100
                    continue
                else:
                    for i in range(len(element1)):
                        replace_char1 = None
                        if lowered_element1[i] != lowered_element2[i]:
                            replace_char1 = lowered_element1[i]
                            replace_char2 = lowered_element2[i]
                            break
                    # Don't replace characters that are numbers or letters, only symbols.
                    # Otherwise 'VITAMIN D3_1' and 'VITAMIN D3_2' get marked as the same.
                    if replace_char1 and \
                       not replace_char1.isalnum() and \
                       not replace_char2.isalnum() and \
                       lowered_element1.replace(replace_char1, replace_char2) == lowered_element2:
                        temp_dict[element2] = 100
                        continue
            
            ratio = fuzzywuzzy.fuzz.ratio(lowered_element1, lowered_element2)
            if ratio >= inclusion_ratio:
                match1 = re.match(r'(.*)_(\d+)', lowered_element1)
                match2 = re.match(r'(.*)_(\d+)', lowered_element2)
                if match1 and match2 and match1.group(2) != match2.group(2):
                    continue
                elif match1:
                    replacement1 = lowered_element2.replace('^', "'")
                    replacement1 = replacement1.replace('_', ',')
                    replacement2 = lowered_element2.replace("'", "^")
                    replacement2 = replacement2.replace(',', '_')
                    replacement3 = lowered_element2.replace("'", "^")
                    replacement3 = replacement3.replace('_', ',')
                    replacement4 = lowered_element2.replace("^", "'")
                    replacement4 = replacement4.replace(',', '_')
                    if lowered_element2 in lowered_element1 or \
                       replacement1 in lowered_element1 or \
                       replacement2 in lowered_element1 or \
                       replacement3 in lowered_element1 or \
                       replacement4 in lowered_element1:
                        continue
                elif match2:
                    replacement1 = lowered_element1.replace('^', "'")
                    replacement1 = replacement1.replace('_', ',')
                    replacement2 = lowered_element1.replace("'", "^")
                    replacement2 = replacement2.replace(',', '_')
                    replacement3 = lowered_element1.replace("'", "^")
                    replacement3 = replacement3.replace('_', ',')
                    replacement4 = lowered_element1.replace("^", "'")
                    replacement4 = replacement4.replace(',', '_')
                    if lowered_element1 in lowered_element2 or \
                       replacement1 in lowered_element2 or \
                       replacement2 in lowered_element2 or \
                       replacement3 in lowered_element2 or \
                       replacement4 in lowered_element2:
                        continue
                
                temp_dict[element2] = ratio
            
            
            
            # ratio = fuzzywuzzy.fuzz.ratio(lowered_element1, element2.lower())
            # if ratio >= 90:
            #     temp_dict[element2] = ratio
            # elif ratio >= 80 and len(element1) == len(element2):
            #     lowered_element2 = element2.lower()
            #     for i in range(len(element1)):
            #         if lowered_element1[i] != lowered_element2[i]:
            #             replace_char1 = lowered_element1[i]
            #             replace_char2 = lowered_element2[i]
            #             break
            #     if lowered_element1.replace(replace_char1, replace_char2) == lowered_element2:
            #         temp_dict[element2] = 100
            #     # replace_chars = [(char,lowered_element2[i]) for i, char in enumerate(lowered_element1) if char not in lowered_element2][0]
            #     # if lowered_element1.replace(replace_chars[0], replace_chars[1]) == lowered_element2:
            #     #     temp_dict[element2] = 100
            #     else:
            #         replacement1 = lowered_element1.replace('^', "'")
            #         replacement1 = replacement1.replace('_', ',')
            #         replacement2 = lowered_element1.replace("'", "^")
            #         replacement2 = replacement2.replace(',', '_')
            #         replacement3 = lowered_element1.replace("'", "^")
            #         replacement3 = replacement3.replace('_', ',')
            #         replacement4 = lowered_element1.replace("^", "'")
            #         replacement4 = replacement4.replace(',', '_')
            #         if replacement1 == lowered_element2 or \
            #             replacement2 == lowered_element2 or \
            #             replacement3 == lowered_element2 or \
            #             replacement4 == lowered_element2:
            #             temp_dict[element2] = 100
        
        if temp_dict:
            fuzz_ratios[element1] = pandas.Series(temp_dict).sort_values()
    return fuzz_ratios


def fuzzy_match(fuzz_ratios):
    """
    """
    matches = {}
    matched = {}
    matches_changed = True
    while matches_changed:
        matches_changed = False
        for key_to_match, ratios in fuzz_ratios.items():
            for i in range(len(ratios)):
                match = ratios.index[i]
                ratio = ratios.iloc[i]
                # Only overwrite the match if the ratio is higher for 
                # the matched key and the matched  
                # key doesn't already have a stronger match.
                if match in matches:
                    if ratio > matches[match]['ratio']:
                        if key_to_match in matched:
                            if matched[key_to_match]['ratio'] < ratio:
                                del matched[matches[match]['match']]
                                matches[match] = {'match':key_to_match, 'ratio':ratio}
                                del matches[matched[key_to_match]['match']]
                                matched[key_to_match] = {'match':match, 'ratio':ratio}
                                matches_changed = True
                                break
                        else:
                            del matched[matches[match]['match']]
                            matches[match] = {'match':key_to_match, 'ratio':ratio}
                            matched[key_to_match] = {'match':match, 'ratio':ratio}
                            matches_changed = True
                            break
                else:
                    if key_to_match in matched:
                        if matched[key_to_match]['ratio'] < ratio:
                            matches[match] = {'match':key_to_match, 'ratio':ratio}
                            del matches[matched[key_to_match]['match']]
                            matched[key_to_match] = {'match':match, 'ratio':ratio}
                            matches_changed = True
                            break
                    else:
                        matches[match] = {'match':key_to_match, 'ratio':ratio}
                        matched[key_to_match] = {'match':match, 'ratio':ratio}
    return matches, matched


list_headers = {"Factors" : {"optional": True, 
                             "all_spellings": ["Factors", "factors"],
                             "index_after_start": 3}, 
                "Samples" : {"optional": False, 
                             "all_spellings": ["Samples", "samples"],
                             "index_after_start": 2}, 
                "metabolite_name" : {"optional": False, 
                                     "all_spellings": ["metabolite_name", "metabolite name", "Metabolite_name", "Metabolite_Name", "Metabolite Name", "Compound", "compound"],
                                     "index_after_start": 2}, 
                "Bin range(ppm)" : {"optional": True, 
                                    "all_spellings": ["Bin range(ppm)", "bin range(ppm)"],
                                    "index_after_start": 2}}

# list headers including misspellings
all_list_headers = [spelling for header_dict in list_headers.values() for spelling in header_dict['all_spellings']]

result_file_keys = ["UNITS", "Has m/z", "Has RT", "RT units"]


def repair(lines, verbose, standardize, last_split, add_data_replace='=', factor_replace=':'):
    """
    """
    # with open(path, encoding="utf-8") as f:
    #     file_text = f.read()
    # file_text = file_text.replace('\n\n', '\n')
    # lines = file_text.split('\n')
    
    
    # # Find METABOLITE_DATA start and end indexes
    # metabolite_data_start_index = [i for i, line in enumerate(lines) if "METABOLITE_DATA_START" in line]
    # if metabolite_data_start_index:
    #     metabolite_data_start_index = metabolite_data_start_index[0]
    
    # metabolite_data_end_index = [i for i, line in enumerate(lines) if "METABOLITE_DATA_END" in line]
    # if metabolite_data_end_index:
    #     metabolite_data_end_index = metabolite_data_end_index[0]
    
    # # Find samples.
    # if metabolite_data_start_index and \
    #    (lines[metabolite_data_start_index+1].startswith("Samples") or \
    #    lines[metabolite_data_start_index+1].startswith("metabolite_name") or \
    #    lines[metabolite_data_start_index+1].startswith("Metabolite Name")):
    #     samples = lines[metabolite_data_start_index+1].split('\t')[1:]
    #     while samples and not samples[-1]:
    #         samples.pop()
    
    
    # # Find METABOLITES start and end indexes
    # metabolites_start_index = [i for i, line in enumerate(lines) if "METABOLITES_START" in line]
    # if metabolites_start_index:
    #     metabolites_start_index = metabolites_start_index[0]
    
    # metabolites_end_index = [i for i, line in enumerate(lines) if "METABOLITES_END" in line]
    # if metabolites_end_index:
    #     metabolites_end_index = metabolites_end_index[0]
    
    # # Find metabolite headers.
    # if metabolites_start_index and \
    #    (lines[metabolites_start_index+1].startswith("Samples") or \
    #    lines[metabolites_start_index+1].startswith("metabolite_name")):
    #     metabolite_headers = lines[metabolites_start_index+1].split('\t')[1:]
    #     while metabolite_headers and not metabolite_headers[-1]:
    #         metabolite_headers.pop()
    
    # # Find EXTENDED_METABOLITE_DATA start and end indexes
    # extended_metabolite_data_start_index = [i for i, line in enumerate(lines) if "EXTENDED" in line and "_METABOLITE_DATA_START" in line]
    # if extended_metabolite_data_start_index:
    #     extended_metabolite_data_start_index = extended_metabolite_data_start_index[0]
    
    # extended_metabolite_data_end_index = [i for i, line in enumerate(lines) if "EXTENDED" in line and "_METABOLITE_DATA_END" in line]
    # if extended_metabolite_data_end_index:
    #     extended_metabolite_data_end_index = extended_metabolite_data_end_index[0]
    
    # # Find extended metabolite headers.
    # if extended_metabolite_data_start_index and \
    #    (lines[extended_metabolite_data_start_index+1].startswith("Samples") or \
    #    lines[extended_metabolite_data_start_index+1].startswith("metabolite_name")):
    #     extended_metabolite_headers = lines[extended_metabolite_data_start_index+1].split('\t')[1:]
    #     while extended_metabolite_headers and not extended_metabolite_headers[-1]:
    #         extended_metabolite_headers.pop()
    
    # # Find BINNED_DATA start and end indexes
    # binned_data_start_index = [i for i, line in enumerate(lines) if "BINNED_DATA_START" in line]
    # if binned_data_start_index:
    #     binned_data_start_index = binned_data_start_index[0]
    
    # binned_data_end_index = [i for i, line in enumerate(lines) if "BINNED_DATA_END" in line]
    # if binned_data_end_index:
    #     binned_data_end_index = binned_data_end_index[0]
    
    # # Find binned headers.
    # if binned_data_start_index and \
    #    lines[binned_data_start_index+1].startswith("Bin range(ppm)"):
    #     binned_headers = lines[binned_data_start_index+1].split('\t')[1:]
    #     while binned_headers and not binned_headers[-1]:
    #         binned_headers.pop()
    
    
    
    
    
    indexes_to_delete = []
    indexes_to_skip = []
    lines_to_add = []
    current_section = None
    current_prefix = ""
    current_section_index = None
    metabolite_data_start_index = None
    metabolites_start_index = None
    extended_metabolite_data_start_index = None
    binned_data_start_index = None
    metabolite_data_end_index = None
    metabolites_end_index = None
    extended_metabolite_data_end_index = None
    binned_data_end_index = None
    ssf_samples = 0
    for i, line in enumerate(lines):
        if i in indexes_to_skip:
            continue
        
        # Upkeep to keep track of where we are in the file.
        if line.startswith('#'):
            current_section = line[1:].strip()
            current_prefix = mwtab.MWTabFile.prefixes.get(current_section, "")
            current_section_index = i
                
        elif "EXTENDED" not in line and "METABOLITE_DATA_START" in line:
            metabolite_data_start_index = i
        elif "METABOLITES_START" in line:
            metabolites_start_index = i
        elif "EXTENDED" in line and "_METABOLITE_DATA_START" in line:
            extended_metabolite_data_start_index = i
        elif "BINNED_DATA_START" in line:
            binned_data_start_index = i
        
        elif "EXTENDED" not in line and "METABOLITE_DATA_END" in line:
            metabolite_data_end_index = i
        elif "METABOLITES_END" in line:
            metabolites_end_index = i
        elif "EXTENDED" in line and "_METABOLITE_DATA_END" in line:
            extended_metabolite_data_end_index = i
        elif "BINNED_DATA_END" in line:
            binned_data_end_index = i
        
        if line.startswith('SUBJECT_SAMPLE_FACTORS'):
            ssf_samples += 1
            
            line_items = line.split("\t")
            
            if line_items[3]:
                factor_pairs = line_items[3].split(" | ")
                factors_dict = {}
                key_count = {}
                for j, pair in enumerate(factor_pairs):
                    factor_data = pair.split(":")
                    if len(factor_data) > 2:
                        if verbose and factor_replace == ':':
                            print("Factor item {} for subject "
                                  "sample factor {} has too many colons, ':'. "
                                  "You can use the \"--last-split\" and \"--factor-replace\" "
                                  "options with this command to fix this issue.".format(j, ssf_samples))
                            print("Factors line:   {}".format(line_items[3]))
                        if last_split:
                            key = factor_replace.join(factor_data[0:-1])
                            value = factor_data[-1]
                        else:
                            key = factor_data[0]
                            value = factor_replace.join(factor_data[1:])
                    else:
                        key = factor_data[0]
                        value = factor_data[1]
                    
                    if key in key_count:
                        factors_dict["{}_{}".format(key, key_count[key])] = value
                        key_count[key] += 1
                    else:
                        key_count[key] = 1
                        factors_dict[key] = value
                factors_data = ["{}:{}".format(key, value) for key, value in factors_dict.items()]
                if any([value > 1 for value in key_count.values()]):
                    print("Duplicate keys in factors.  Line number: " + str(i))
                line_items[3] = " | ".join(factors_data)
                
            
            if line_items[4]:
                # Some additional data items are separated by 2 spaces instead of '; '
                if any(len(factor_item.split('=')) > 2 for factor_item in line_items[4].split('; ')):
                    line_items[4] = line_items[4].replace('  ', '; ')
                
                # Fix duplicate keys and replace '=' if more than 1.
                additional_data_dict = {}
                key_count = {}
                for j, factor_item in enumerate(line_items[4].split("; ")):
                    add_data = factor_item.split("=")
                    if len(add_data) > 2:
                        if verbose and add_data_replace == '=':
                            print("Additional sample data item {} for subject "
                                  "sample factor {} has too many equal signs, '='. "
                                  "You can use the \"--last-split\" and \"--add-replace\" "
                                  "options with this command to fix this issue.".format(j, ssf_samples))
                            print("Additional data line:   {}".format(line_items[4]))
                        if last_split:
                            key = add_data_replace.join(add_data[0:-1])
                            value = add_data[-1]
                        else:
                            key = add_data[0]
                            value = add_data_replace.join(add_data[1:])
                    else:
                        key = add_data[0]
                        value = add_data[1]
                    
                    if key in key_count:
                        additional_data_dict["{}_{}".format(key, key_count[key])] = value
                        key_count[key] += 1
                    else:
                        key_count[key] = 1
                        additional_data_dict[key] = value
                additional_sample_data = ["{}={}".format(key, value) for key, value in additional_data_dict.items()]
                if any([value > 1 for value in key_count.values()]):
                    print("Duplicate keys in additional data.  Line number: " + str(i))
                line_items[4] = "; ".join(additional_sample_data)
                
            lines[i] = '\t'.join(line_items)
            continue
        
        
        # Remove tabs from the workbench header.
        if line.startswith("#METABOLOMICS WORKBENCH") and '\t' in line:
            lines[i] = re.sub(r'\t+', ' ', line)
            print("Tabs removed from WORKBENCH header.  Line number: " + str(i))
            continue
        
        # Fix it when the header accidently is on 2 different lines.
        if i > 0 and lines[i-1].startswith("#METABOLOMICS WORKBENCH") and '\t' not in line:
                lines[i-1] = lines[i-1] + ' ' + line.lstrip()
                indexes_to_delete.append(i)
                print("WORKBENCH header on multiple lines.  Line number: " + str(i))
                continue
        
        # If VERSION doesn't have a value then give it one.
        if match_re := re.match(r"VERSION(.*)", line):
            version_value = match_re.group(1).strip()
            if version_value == '':
                lines[i] = "VERSION             \t1"
                print("Missing version number.  Line number: " + str(i))
                continue
            else:
                lines[i] = "VERSION             \t" + version_value
                continue
            # if '\t' not in match_re.group(1):
            #     lines[i] = "VERSION             \t" + match_re.group(1).lstrip()
            #     continue
        
        # Fix tabs in the sub_section names.
        if (match_re := re.match(r'(\s*[A-Z_]+\s*[A-Z_]*:?\s*[A-Z_]+)( +\t.*)', line)) and '\t' in match_re.group(1):
            lines[i] = "".join(match_re.group(1).split()) + match_re.group(2)
            print("Tab(s) in sub section name.  Line number: " + str(i))
            continue
        
        # Fix tabs in #SECTIONS.
        if line.startswith('#') and \
           '\t' in line and \
           "RESULTS_FILE" not in line and\
           "SUBJECT_SAMPLE_FACTORS" not in line:
            lines[i] = "".join(line.split())
            print("Tab(s) in section title, (#SECTION line).  Line number: " + str(i))
            continue
        
        # Fix names for list header lines.
        if metabolite_data_start_index:
            if i - metabolite_data_start_index == 1:
                # This can't be done for samples because they have to match up to the samples in SSF.
                # samples = line.split('\t')[1:]
                # while samples and not samples[-1]:
                #     samples.pop()
                
                # header_count = {}
                # for j, header in enumerate(samples):
                #     if header in header_count:
                #         samples[j] = f"{header}_{header_count[header]}"
                #         header_count[header] += 1
                #     else:
                #         header_count[header] = 1
                
                # lines[i] = '\t'.join(["Samples"] + samples)
                # continue
                
                # num_of_metabolite_data_headers = len(line.split('\t'))
                if not line.startswith("Samples"):
                    if '\t' in line:
                        lines[i] = re.sub(r"[^\t]*\t", "Samples\t", line, count=1)
                    else:
                        lines[i] = "Samples"
                    print("METABOLITE_DATA missing 'Samples' in header.  Line number: " + str(i))
                    continue
                
                if line.startswith("Samples\tSamples"):
                    lines[i] = line.replace("Samples\tSamples", "Samples")
                    print("METABOLITE_DATA has too many 'Samples' in header.  Line number: " + str(i))
                continue
            
            elif i - metabolite_data_start_index == 2:
                needs_continue = False
                if not line.startswith("Factors") and not line_is_numbers(line):
                    lines[i] = re.sub(r"[^\t]*\t", "Factors\t", line, count=1)
                    line = lines[i]
                    needs_continue = True
                    print("METABOLITE_DATA missing 'Factors' in factor line.  Line number: " + str(i))
                
                if line.startswith("Factors"):
                    new_factors = []
                    sample_factors = line.split('\t')[1:]
                    while sample_factors and not sample_factors[-1]:
                        sample_factors.pop()
                    
                    for k, sample_factor_sets in enumerate(sample_factors):
                        factor_pairs = sample_factor_sets.split(" | ")
                        factors_dict = {}
                        key_count = {}
                        for j, pair in enumerate(factor_pairs):
                            factor_data = pair.split(":")
                            if len(factor_data) > 2:
                                if verbose and factor_replace == ':':
                                    print("The Factor line in the DATA section has a sample, #{}, with a factor item, #{}, "
                                          "that has too many colons, ':'. "
                                          "You can use the \"--last-split\" and \"--factor-replace\" "
                                          "options with this command to fix this issue.".format(k, j))
                                    print("Factor item:   {}".format(pair))
                                if last_split:
                                    key = factor_replace.join(factor_data[0:-1])
                                    value = factor_data[-1]
                                else:
                                    key = factor_data[0]
                                    value = factor_replace.join(factor_data[1:])
                            else:
                                key = factor_data[0]
                                value = factor_data[1]
                            
                            if key in key_count:
                                factors_dict["{}_{}".format(key, key_count[key])] = value
                                key_count[key] += 1
                            else:
                                key_count[key] = 1
                                factors_dict[key] = value
                        factors_data = ["{}:{}".format(key, value) for key, value in factors_dict.items()]
                        if any([value > 1 for value in key_count.values()]):
                            print("Duplicate keys in DATA factors.  Line number: " + str(i))
                        new_factors.append(" | ".join(factors_data))
                    lines[i] = '\t'.join(["Factors"] + new_factors)
                    needs_continue = True
                
                if needs_continue:
                    continue
        
        if metabolites_start_index:
            if i - metabolites_start_index == 1:
                metabolite_headers = line.split('\t')
                # Rare case where the header line is missing. Only seen it happen once and with only 1 column.
                if len(metabolite_headers) == 1 and not any(metabolite_headers[0] == header for header in all_list_headers):
                    lines_to_add.append(['metabolite_name', 'METABOLITES'])
                    print("METABOLITES missing header line.  Line number: " + str(i))
                    continue
                
                metabolite_headers = metabolite_headers[1:]
                # num_of_metabolites_headers = len(metabolite_headers) + 1
                while metabolite_headers and not metabolite_headers[-1]:
                    metabolite_headers.pop()
                
                if "Metabolite" in metabolite_headers:
                    metabolite_headers = [header if header != "Metabolite" else "Metabolite Name" for header in metabolite_headers]
                    print("Extra 'Metabolite' header in METABOLITES.  Line number: " + str(i))
                
                header_count = {}
                for j, header in enumerate(metabolite_headers):
                    if standardize:
                        # Check if headers are recognized variations and change it to the standardized name.
                        if not any(standardized_header == header for standardized_header in METABOLITES_REGEXS.keys()):
                            for standardized_header in METABOLITES_REGEXS.keys():
                                if any(re.match(p, header) for p in METABOLITES_REGEXS[standardized_header]):
                                    header = standardized_header
                                    print("METABOLITES has non standard header(s).  Line number: " + str(i))
                                    # TODO fix this. The header doesn't actually change as is. AN000063
                    
                    if header in header_count:
                        metabolite_headers[j] = "{}_{}".format(header, header_count[header])
                        header_count[header] += 1
                    else:
                        header_count[header] = 1
                
                if any([value > 1 for value in header_count.values()]):
                    print("Duplicate headers in METABOLITES.  Line number: " + str(i))
                lines[i] = '\t'.join(["metabolite_name"] + metabolite_headers)
                continue
                
                # if not line.startswith("metabolite_name"):
                #     lines[i] = re.sub(r"[^\t]*\t", "metabolite_name\t", line, count=1)
                #     continue
        
        if extended_metabolite_data_start_index:
            if i - extended_metabolite_data_start_index == 1:
                extended_headers = line.split('\t')[1:]
                # num_of_extended_headers = len(extended_headers) + 1
                while extended_headers and not extended_headers[-1]:
                    extended_headers.pop()
                
                if "Metabolite" in extended_headers:
                    extended_headers = [header if header != "Metabolite" else "Metabolite Name" for header in extended_headers]
                    print("Extra 'Metabolite' header in EXTENDED.  Line number: " + str(i))
                
                header_count = {}
                for j, header in enumerate(extended_headers):
                    if header in header_count:
                        extended_headers[j] = "{}_{}".format(header, header_count[header])
                        header_count[header] += 1
                    else:
                        header_count[header] = 1
                
                if any([value > 1 for value in header_count.values()]):
                    print("Duplicate headers in EXTENDED.  Line number: " + str(i))
                lines[i] = '\t'.join(["metabolite_name"] + extended_headers)
                continue
                
                # if not line.startswith("metabolite_name"):
                #     lines[i] = re.sub(r"[^\t]*\t", "metabolite_name\t", line, count=1)
                #     continue
        
        if binned_data_start_index:
            if i - binned_data_start_index == 1:
                # This can't be done for binned_headers because they have to match up to the samples in SSF.
                # binned_headers = line.split('\t')[1:]
                # while binned_headers and not binned_headers[-1]:
                #     binned_headers.pop()
                
                # header_count = {}
                # for j, header in enumerate(binned_headers):
                #     if header in header_count:
                #         binned_headers[j] = f"{header}_{header_count[header]}"
                #         header_count[header] += 1
                #     else:
                #         header_count[header] = 1
                
                # lines[i] = '\t'.join(["Bin range(ppm)"] + binned_headers)
                # continue
                
                # num_of_binned_headers = len(line.split('\t'))
                if not line.startswith("Bin range(ppm)"):
                    lines[i] = re.sub(r"[^\t]*\t", "Bin range(ppm)\t", line, count=1)
                    print("BINNNED DATA missing 'Bin range(ppm) header.  Line number: " + str(i))
                    continue
                continue
        
        # Find header lines in the wrong place and delete them.
        first_tab = line.split('\t')[0]
        if any(first_tab == header for header in all_list_headers):
        # if any(line.startswith(header) for header in all_list_headers):
            # If it's the wrong index it should have taken one of the if's above this.
            # if i - metabolite_data_start_index == 1 or\
            #    i - metabolite_data_start_index == 2 or\
            #    i - metabolites_start_index == 1 or\
            #    i - extended_metabolite_data_start_index == 1 or\
            #    i - binned_data_start_index == 1:
            #     continue
            # We assume that if this line is not where we expect it to be then it shouldn't be there.
            indexes_to_delete.append(i)
            print("Header line in the wrong place.  Line number: " + str(i))
        
        # Fix the issue where tabs inside of double quotes messes things up in data blocks.
        if data_block_element_fix(metabolite_data_start_index, metabolite_data_end_index, line, i, lines):
            print("Tabs inside double quotes (METABOLITE DATA).  Line number: " + str(i))
            continue
        if data_block_element_fix(metabolites_start_index, metabolites_end_index, line, i, lines):
            print("Tabs inside double quotes (METABOLITES).  Line number: " + str(i))
            continue
        if data_block_element_fix(extended_metabolite_data_start_index, extended_metabolite_data_end_index, line, i, lines):
            print("Tabs inside double quotes (EXTENDED).  Line number: " + str(i))
            continue
        if data_block_element_fix(binned_data_start_index, binned_data_end_index, line, i, lines):
            print("Tabs inside double quotes (BINNED DATA).  Line number: " + str(i))
            continue
        
        # Fix bad prefix.
        if not line.startswith(current_prefix) and not line.startswith('#'):
            # The line could have a bad prefix or no prefix.
            # Assume that if a ':' is found close enough to the start of the line it is a bad prefix.
            # "_1 CH: ..." has been seen before, so using index 6 as the "close enough".
            if ':' in line and line.index(':') < 6:
                lines[i] = re.sub(r"[^:]*:", current_prefix, line, count=1)
                print("Bad prefix.  Line number: " + str(i))
            # Assume if there is no tab that the line should be added to the one above it.
            elif '\t' not in line:
                lines[i-1] = lines[i-1] + line.lstrip()
                indexes_to_delete.append(i)
                print("Unexpected newline.  Line number: " + str(i)) # Line moved up.
            else:
                lines[i] = current_prefix + line
                print("Missing prefix.  Line number: " + str(i))
            continue
        
        # Remove #FACTORS lines.
        if line.strip() == "#FACTORS":
            indexes_to_delete.append(i)
            print("Invalid '#FACTORS' section.  Line number: " + str(i))
            continue
        
        # Remove extra #END lines.
        if line.strip() == "#END" and any(["#END" in prev_line for prev_line in lines[i-5:i]]):
            indexes_to_delete.append(i)
            print("Extra '#END' line.  Line number: " + str(i))
            continue
        
        # Fix various problems with RESULTS_FILE lines.
        if "RESULTS_FILE" in line:
            problem_found = False
            
            if line.startswith('#'):
                line = line[1:]
                problem_found = True
            
            if '\t\t' in line:
                line = line.replace('\t\t', '\t')
                problem_found = True
            
            if "NM_RESULTS_FILE" in line:
                line = line.replace("NM_RESULTS_FILE", "NMR_RESULTS_FILE")
                problem_found = True
            
            if "NMR_RESULTS_FILE" in line:
                if "NM:" not in line:
                    line = "NM:" + line
                    problem_found = True
                if "NM:" in line and line[0:3] != "NM:":
                    line = re.sub(r'[^:]*:', 'NM:', line, 1)
                    problem_found = True
                if current_prefix != 'NM:':
                    lines_to_add.append([line, 'NM'])
                    indexes_to_delete.append(i)
                    problem_found = True
            
            if "MS_RESULTS_FILE" in line:
                if "MS:" not in line:
                    line = "MS:" + line
                    problem_found = True
                if "MS:" in line and line[0:3] != "MS:":
                    line = re.sub(r'[^:]*:', 'MS:', line, 1)
                    problem_found = True
                if current_prefix != 'MS:':
                    lines_to_add.append([line, 'MS'])
                    indexes_to_delete.append(i)
                    problem_found = True
                  
            for result_key in result_file_keys:
                if result_key in line and not re.match(r'.*\w+ *\t *{}.*'.format(result_key), line):
                    line = line.replace(result_key, '\t' + result_key)
                    problem_found = True
                
            lines[i] = line
            if problem_found:
                print("Malformed RESULTS_FILE.  Line number: " + str(i))
            continue
        
        # Remove duplicate sub sections and other stuff.
        if match_re := re.match(r'((\w\w):([A-Z_]+))(.*)', line):
            # Look for duplicate sub sections.
            sub_section_indexes = [index + current_section_index for index, line in enumerate(lines[current_section_index:]) if line.startswith(match_re.group(1))]
            duplicate_line_indexes = [index + current_section_index for index, line2 in enumerate(lines[current_section_index:]) if line2 == line]
            # I can't remember why I put the != condition in, so I am commenting it out for now in case I find the reason again later.
            # if len(duplicate_line_indexes) > 1 and len(sub_section_indexes) != len(duplicate_line_indexes):
            if len(duplicate_line_indexes) > 1:
                indexes_to_delete.extend(duplicate_line_indexes[1:])
                indexes_to_skip.extend(duplicate_line_indexes)
                print("Duplicate sub sections.  Line number: " + str(i))
                continue
            
            # Add tab in if not there. 
            if '\t' not in match_re.group(4):
                # Adding the correct spacing isn't necessary because 
                # using the MWTabFile class at the end should fix that, but I am doing it anyway.
                lines[i] = match_re.group(1) + (30 - len(match_re.group(1))) * ' ' + '\t' + match_re.group(4).lstrip()
                print("Incorrect sub section spacing.  Line number: " + str(i))
                continue
        
        # Remove out of place sub section.
        if line.startswith("PR:INSTITUTE") and "ST:INSTITUTE" in line:
            lines[i] = line[0:line.index("ST:INSTITUTE")]
            print("Out of place ST:INSTITUTE.  Line number: " + str(i))
            continue
    
    for i, index in enumerate(sorted(indexes_to_delete)):
        del lines[index - i]
    
    sections_with_start = ['METABOLITES', 'NMR_METABOLITE_DATA', 'MS_METABOLITE_DATA', 'NMR_BINNED_DATA', 'EXTENDED_METABOLITE_DATA']
    for line, section_name in lines_to_add:
        try:
            section_index = [index for index, line2 in enumerate(lines) if line2 == "#" + section_name][0]
        except Exception:
            continue
        if section_name in sections_with_start:
            lines.insert(section_index + 2, line)
        else:
            lines.insert(section_index + 1, line)


    # Try to use the MWTabFile class to get an internal representation and save it out because this will standardize spacing issues.
    # Also fixes broken subsections.
    # If that can't be done then just save out the lines as is.
    file_str = '\n'.join(lines)
    
    tabfile = mwtab.MWTabFile("", compatability_mode=True)
    conversion_successful = False
    try:
        tabfile._build_mwtabfile(file_str)
        conversion_successful = True
    except Exception as e:
        if verbose:
            traceback.print_exception(e, file=sys.stdout)
            print()
        return file_str
    
    if conversion_successful:
        
        data_section_key = list(set(tabfile.keys()) &
                                {"MS_METABOLITE_DATA", "NMR_METABOLITE_DATA", "NMR_BINNED_DATA"})
        if data_section_key:
            data_section_key = data_section_key[0]
            data_metabolites = None
            metabolites_metabolites = None
            tables = {}
            for table_name in ['Data', 'Metabolites', 'Extended']:
                if table_name in tabfile[data_section_key]:
                    temp_list = [duplicates_dict._JSON_DUPLICATE_KEYS__Jobj for duplicates_dict in tabfile[data_section_key][table_name]]
                    data_df = pandas.DataFrame.from_records(temp_list).astype(str)
                    # Fix any html characters.
                    data_df = data_df.map(html.unescape)
                    
                    data_df = clean_columns(data_df)
                    
                    # TODO after adding the code to shift bad values, look for "Unnamed" columns and remove them or at least rename them back to empty string.
                    # Can repeat what is done in the clean_columns function.
                    # The column names should not change after repair.
                    
                    # Look for completely null columns and remove if found.
                    no_name_columns = [column for column in data_df.columns if column == '' or re.match(r'\{\{\{_\d{1,}_\}\}\}', column)]
                    if no_name_columns:
                        # Drop null columns, non-numeric columns, and columns that are the sum of the other columns.
                        no_name_df = data_df.loc[:, no_name_columns]
                        null_columns = no_name_df.isna().all() | (no_name_df == "").all()
                        
                        if table_name == 'Data':
                            non_numeric_columns = ~(no_name_df.apply(is_column_numeric) | null_columns)
                            # If a column is the sum of other columns then drop it.
                            numeric_df = _create_numeric_df(data_df)
                            sum_df = numeric_df.loc[:, [column for column in numeric_df.columns if column not in no_name_columns]].sum(axis=1)
                            numeric_no_name_df = no_name_df.apply(pandas.to_numeric, errors='coerce')
                            # Equal comparison doesn't work because of floating point math.
                            # Calculate percent difference and if it is low enough say they are equal.
                            percent_diff = numeric_no_name_df.sub(sum_df, axis='index').div(sum_df, axis='index').abs()
                            columns_to_drop = null_columns | non_numeric_columns | (percent_diff < 1.0e-6).all()
                        else:
                            columns_to_drop = null_columns
                        
                        columns_to_drop = columns_to_drop[columns_to_drop].index
                        
                        # TODO make sure all of these meet the criteria to be dropped.
                        if len(columns_to_drop) > 0:
                            print("Null column deleted in '" + table_name)
                            data_df.drop(columns = columns_to_drop, inplace=True)
                            if table_name == 'Data':
                                tabfile._samples = [column for column in data_df.columns if column != "Metabolite"]
                            elif table_name == 'Metabolites':
                                tabfile._metabolite_header = [column for column in data_df.columns if column != "Metabolite"]
                            elif table_name == 'Extended':
                                tabfile._extended_metabolite_header = [column for column in data_df.columns if column != "Metabolite"]
                                        
                    # Consolidate values by column.
                    if table_name == 'Metabolites':
                        consolidate_values_to_column(data_df, 'Metabolites', ['kegg', 'id'], r'[cC]\d{5}')
                        consolidate_values_to_column(data_df, 'Metabolites', ['inchi', 'key'], r'[a-zA-Z]{14}-[a-zA-Z]{10}-[a-zA-Z]')
                    
                    # TODO look to see if the metabolites column is duplicated and 
                    # that there is a column name of empty string at the end and it has values, then move all columns left one.
                    # Example: AN004528
                    # Fix missing tabs.
                    if table_name in ['Extended', 'Metabolites']:
                        fix_missing_tab_in_table(data_df, table_name)
                                        
                    # Fill any na values with the empty string.
                    data_df = data_df.fillna('')
                    
                    # Get lists of metabolites for later.
                    if table_name == 'Data' and "Metabolite" in data_df.columns:
                        data_metabolites = list(data_df.loc[:, "Metabolite"])
                    if table_name == 'Metabolites' and "Metabolite" in data_df.columns:
                        metabolites_metabolites = list(data_df.loc[:, "Metabolite"])
                    
                    # TODO fix this so the column names don't have the {{{_\d_}}} in them. I think you can just rename the columns.
                    tabfile[data_section_key][table_name] = [DuplicatesDict(data_dict) for data_dict in data_df.astype(str).to_dict(orient='records')]
                    tables[table_name] = data_df
            
            # Remove duplicate rows from Data.
            # TODO change this so that we look for metabolites in Data that don't match in Metabolites first, then we see if they have a 
            # duplicate and if they do we can drop it from Data. We can maybe drop the fuzzy ratio threshold as well since the numeric 
            # matching gives a stronger chance that it is the same metabolite.
            # data_df = tables['Data']
            # if len(data_df) > 0:
            #     data_df_init_len = len(data_df)
            #     # Have a slightly different criteria for removing some columns from determining duplicates.
            #     no_name_columns = [column for column in data_df.columns if column == '' or re.match(r'\{\{\{_\d{1,}_\}\}\}', column)]
            #     if no_name_columns:
            #         no_name_df = data_df.loc[:, no_name_columns]
            #         null_columns = no_name_df.isna().all() | (no_name_df == "").all() | (no_name_df.isna() | (no_name_df == '0')).all()
            #         columns_to_drop = null_columns[null_columns].index
            #         if len(columns_to_drop) > 0:
            #             dup_df = data_df.drop(columns = columns_to_drop)
            #         else:
            #             dup_df = data_df
            #     else:
            #         dup_df = data_df
                
            #     # Convert every column except the metabolite column to numbers.
            #     numeric_df = dup_df.iloc[:, 1:]
            #     numeric_df = numeric_df.apply(pandas.to_numeric, errors='coerce')
            #     # Add metabolites back and drop duplicates.
            #     numeric_df.loc[:, 'Metabolites'] = dup_df.loc[:, 'Metabolite']
            #     numeric_df = numeric_df.drop_duplicates()
            #     # Remove Metabolites and find duplicate rows.
            #     numeric_df = numeric_df.drop('Metabolites', axis=1)
            #     duplicates_bool = numeric_df.duplicated(keep=False)
            #     duplicates = numeric_df[duplicates_bool]
            #     # Remove rows that are mono value. For example, if a row has all nan's that should not be considered as a duplicate.
            #     duplicates = duplicates[~duplicates.nunique(axis = 1, dropna=False).eq(1)]
            #     if len(duplicates) != 0:
            #         groups = duplicates.groupby(duplicates.columns.tolist(), dropna=False)
            #         data_indexes_to_delete = []
            #         met_indexes_to_delete = []
            #         for name, group in groups:
            #             names = list(data_df.loc[group.index, :].iloc[:, 0])
            #             fuzz_ratios = compute_fuzz_ratios(names, names)
            #             fuzz_ratios = {metabolite:series.drop(metabolite, errors='ignore') for metabolite, series in fuzz_ratios.items()}
            #             met_to_index = {metabolite:group.iloc[i].name for i, (metabolite, series) in enumerate(fuzz_ratios.items()) if len(series) > 0}
            #             # Add all indexes to list to delete except for 1.
            #             if metabolites_metabolites:
            #                 duplicates_in_metabolites = {metabolite:index for metabolite, index in met_to_index.items() if metabolite in metabolites_metabolites}
            #                 duplicates_not_in_metabolites = {metabolite:index for metabolite, index in met_to_index.items() if metabolite not in metabolites_metabolites}
            #                 # See if the duplicates in Data have duplicate metadata in Metabolites.
            #                 met_dups = tables['Metabolites'][tables['Metabolites'].loc[:, 'Metabolite'].isin(duplicates_in_metabolites.keys())]
            #                 met_dups_data = met_dups.drop('Metabolite', axis=1)
            #                 met_dups = met_dups[met_dups_data.duplicated(keep=False)].loc[:, 'Metabolite']
            #                 matching_dups = {metabolite:index for metabolite, index in duplicates_in_metabolites.items() if metabolite in met_dups}
            #                 # non_matching_dups = {metabolite:index for metabolite, index in duplicates_in_metabolites if metabolite not in met_dups}
            #                 # Only remove duplicates that are not in metabolites or are also duplicates in metabolites.
            #                 duplicate_indexes = list(duplicates_not_in_metabolites.values()) + \
            #                                     list(matching_dups.values())
            #                 met_duplicate_indexes = [met_dups[metabolite].index[0] for metabolite, index in matching_dups.items() if metabolite in met_dups]
            #                 if len(duplicate_indexes) == len(met_to_index):
            #                     duplicate_indexes = duplicate_indexes[:-1]     
            #                     met_duplicate_indexes = met_duplicate_indexes[:-1]
            #             else:
            #                 duplicate_indexes = list(met_to_index.values())[:-1]
            #             data_indexes_to_delete += duplicate_indexes
            #             met_indexes_to_delete += met_duplicate_indexes
            #         numeric_df = numeric_df[~numeric_df.index.isin(data_indexes_to_delete)]
            #         data_df = data_df.loc[numeric_df.index, :]
            #         met_df = tables['Metabolites'].loc[[index for index in tables['Metabolites'].index if index not in met_indexes_to_delete], :]
            #         if data_df_init_len > len(data_df):
            #             print("Duplicate rows removed from DATA section.")
            #             tabfile[data_section_key]['Data'] = [DuplicatesDict(data_dict) for data_dict in data_df.astype(str).to_dict(orient='records')]
            #             tabfile[data_section_key]['Metabolites'] = [DuplicatesDict(data_dict) for data_dict in met_df.astype(str).to_dict(orient='records')]
            
            # Match up metabolites in DATA section and METABOLITES sections for metabolites with slightly different names.
            # For example, 2^,4^-compound and 2'_4'-compound are the same, but for some reason one name is in DATA and another in METABOLITES.
            if data_metabolites and metabolites_metabolites:
                in_data_not_met = [metabolite for metabolite in data_metabolites if metabolite not in metabolites_metabolites]
                in_met_not_data = [metabolite for metabolite in metabolites_metabolites if metabolite not in data_metabolites]
                data_df = tables['Data']
                
                if in_data_not_met and in_met_not_data:
                    data_fuzz_ratios = compute_fuzz_ratios(in_data_not_met, in_met_not_data)
                    met_fuzz_ratios = compute_fuzz_ratios(in_met_not_data, in_data_not_met)
                    
                    # data_matches has keys that are the metabolites in in_met_not_data
                    # matched_data has keys that are the metabolites in in_data_not_met
                    data_matches, matched_data = fuzzy_match(data_fuzz_ratios)
                    met_matches, matched_met = fuzzy_match(met_fuzz_ratios)
                    
                    data_pairs = {(metabolite, match['match']) for metabolite, match in data_matches.items()}
                    met_pairs = {(match['match'], metabolite) for metabolite, match in met_matches.items()}
                    
                    # All pairs that are the same are good to go.
                    definite_matches = data_pairs.intersection(met_pairs)
                    # For pairs that aren't shared between the 2 sets if both pair 
                    # names are unique to that set of pairs then it can be added 
                    # because there is no conflict with the other set of pairs.
                    data_pairs_not_in_met = met_pairs - data_pairs
                    conflicted_data_pairs_metabolites = {metabolite for pair in data_pairs_not_in_met for metabolite in pair}
                    met_pairs_not_in_data = data_pairs - met_pairs
                    conflicted_met_pairs_metabolites = {metabolite for pair in met_pairs_not_in_data for metabolite in pair}
                    for data_pair in data_pairs_not_in_met:
                        if data_pair[0] not in conflicted_met_pairs_metabolites \
                           and data_pair[1] not in conflicted_met_pairs_metabolites:
                            definite_matches.add(data_pair)
                    for met_pair in met_pairs_not_in_data:                            
                        if met_pair[0] not in conflicted_data_pairs_metabolites \
                           and met_pair[1] not in conflicted_data_pairs_metabolites:
                            definite_matches.add(met_pair)
                    
                    definite_matches = {pair[1]:pair[0] for pair in definite_matches}
                    if definite_matches:
                        # temp_list = [duplicates_dict._JSON_DUPLICATE_KEYS__Jobj for duplicates_dict in tabfile[data_section_key]['Data']]
                        # data_df = pandas.DataFrame.from_records(temp_list)
                        ## data_df = pandas.DataFrame.from_records(tabfile[data_section_key]['Data'])
                        data_df.loc[:, "Metabolite"] = data_df.loc[:, "Metabolite"].replace(definite_matches)
                        tabfile[data_section_key]['Data'] = [DuplicatesDict(data_dict) for data_dict in data_df.astype(str).to_dict(orient='records')]
                        print("The following metabolites in the 'DATA' section were changed "
                              "to the indicated metabolite name that was fuzzy matched from "
                              "the 'METABOLITES' section:")
                        try:
                            print("\n".join([key + "   ->   " + value for key, value in definite_matches.items()]))
                        except UnicodeEncodeError:
                            print("\n".join([key + "   ->   " + value for key, value in definite_matches.items()]).encode('utf-8'))
                    
                    
                    # For left over Data metabolties not in Metabolites see if they have a duplicate that is in Metabolites and delete if so.
                    leftover_data_mets = [metabolite for metabolite in in_data_not_met if metabolite not in definite_matches]
                    if leftover_data_mets:
                        # dup_df = data_df[~data_df.loc[:, 'Metabolite'].isin(leftover_data_mets)]
                        
                        # Have a slightly different criteria for removing some columns from determining duplicates.
                        no_name_columns = [column for column in data_df.columns if column == '' or re.match(r'\{\{\{_\d{1,}_\}\}\}', column)]
                        if no_name_columns:
                            no_name_df = data_df.loc[:, no_name_columns]
                            null_columns = no_name_df.isna().all() | (no_name_df == "").all() | (no_name_df.isna() | (no_name_df == '0')).all()
                            columns_to_drop = null_columns[null_columns].index
                            if len(columns_to_drop) > 0:
                                dup_df = data_df.drop(columns = columns_to_drop)
                            else:
                                dup_df = data_df
                        else:
                            dup_df = data_df
                        
                        # Convert every column except the metabolite column to numbers.
                        numeric_df = dup_df.iloc[:, 1:]
                        numeric_df = numeric_df.apply(pandas.to_numeric, errors='coerce')
                        # Add metabolites back and drop duplicates.
                        numeric_df.loc[:, 'Metabolites'] = dup_df.loc[:, 'Metabolite']
                        numeric_df = numeric_df.drop_duplicates()
                        # Remove Metabolites and find duplicate rows.
                        numeric_df = numeric_df.drop('Metabolites', axis=1)
                        duplicates_bool = numeric_df.duplicated(keep=False)
                        duplicates = numeric_df[duplicates_bool]
                        # Remove rows that are mono value. For example, if a row has all nan's that should not be considered as a duplicate.
                        duplicates = duplicates[~duplicates.nunique(axis = 1, dropna=False).eq(1)]
                        
                        # Identify families of metabolites.
                        # A family is something like 'VITAMIN B12' and 'VITAMIN B12_1'
                        # data_metabolites could be changed due to renaming above, so rebuild it.
                        data_metabolites = list(data_df.loc[:, "Metabolite"])
                        roots = [met for met in data_metabolites if not re.match(r'.*_\d+$', met)]
                        root_to_mets = {root:[root] for root in roots}
                        met_to_root = {root:root for root in roots}
                        for met in data_metabolites:
                            if met in roots:
                                continue
                            for root in roots:
                                if root in met:
                                    met_to_root[met] = root
                                    root_to_mets[root].append(met)
                                    
                        # Remove duplicates that are family members.
                        duplicate_mets = data_df.loc[duplicates.index, 'Metabolite']
                        dup_indexes_to_keep = []
                        if len(duplicates) != 0:
                            groups = duplicates.groupby(duplicates.columns.tolist(), dropna=False)
                            for name, group in groups:
                                names = list(data_df.loc[group.index, :].iloc[:, 0])
                                roots_to_names = {}
                                for name2 in names:
                                    root = met_to_root[name2]
                                    if root in roots_to_names:
                                        roots_to_names[root].append(name2)
                                    else:
                                        roots_to_names[root] = [name2]
                                names_to_keep = [names[0] for root, names in roots_to_names.items() if len(names) == 1]
                                dup_indexes_to_keep += list(data_df[data_df.loc[:, 'Metabolite'].isin(names_to_keep)].index)
                        duplicate_mets = duplicate_mets.loc[dup_indexes_to_keep]
                        
                        # Remove leftover metabolites in the data section that aren't in the metabolites section and have a duplicate.
                        leftover_mets_in_dups = [metabolite for metabolite in leftover_data_mets if metabolite in duplicate_mets.values]
                        if leftover_mets_in_dups:
                            data_df = data_df[~data_df.loc[:, 'Metabolite'].isin(leftover_mets_in_dups)]
                            print("Duplicate rows removed from DATA section.")
                            tabfile[data_section_key]['Data'] = [DuplicatesDict(data_dict) for data_dict in data_df.astype(str).to_dict(orient='records')]
                        
                                
        return tabfile.writestr('mwtab')
    else:
        return file_str
    # try:
    #     tabfile._build_mwtabfile(file_str)
    #     return tabfile.writestr('mwtab')
    # except Exception as e:
    #     if verbose:
    #         traceback.print_exception(e, file=sys.stdout)
    #         print()
    #     return file_str




def clean_columns(data_df):
    """
    """
    for i, column in enumerate([column for column in data_df.columns if column != 'Metabolite']):
        data_df.loc[:, column] = data_df.loc[:, column].str.strip()
        data_df.loc[:, column] = data_df.loc[:, column].str.strip('\u200e')
        # Replace wierd unicode and HTML characters.
        data_df.loc[:, column] = data_df.loc[:, column].str.replace('\u3000', ' ')
        data_df.loc[:, column] = data_df.loc[:, column].str.replace('&gt;', '>')
        data_df.loc[:, column] = data_df.loc[:, column].str.replace('&lt;', '<')
        data_df.loc[:, column] = data_df.loc[:, column].str.replace('âˆ’', '-')
        data_df.loc[:, column] = data_df.loc[:, column].str.replace('Â ', '')
        data_df.loc[:, column] = data_df.loc[:, column].str.replace('−', '-')
        data_df.loc[:, column] = data_df.loc[:, column].str.replace('–', '-')
        # Replace NA values.
        # TODO look to see if 0's are being used as NA values and replace them.
        lowered_column = column.lower()
        if lowered_column.endswith(' id') or lowered_column.endswith('_id') or 'inchi' in lowered_column or \
           'kegg' in lowered_column or 'hmdb' in lowered_column:
            condition = data_df.loc[:, column] == '0'
            data_df.loc[condition, column] = ''
        na_values = data_df.loc[:, column].isin(NA_VALUES)
        data_df.loc[na_values, column] = ''
        # Looking for numeric columns that have commas instead of decimals and replacing the commas with decimals.
        temp_column = data_df.loc[:, column]
        comma_count = temp_column.str.count(',')
        single_commas = comma_count == 1
        no_spaces_or_periods = ~temp_column.str.contains('\.| ', regex=True, na=False)
        replace_candidates = (single_commas & no_spaces_or_periods) | (temp_column == '')
        if replace_candidates.all():
            temp_column = temp_column[replace_candidates].str.replace(',', '.')
            to_numeric = pandas.to_numeric(temp_column, errors='coerce')
            commas_to_replace = replace_candidates & ~to_numeric.isna()
            data_df.loc[:, column] = data_df.loc[commas_to_replace, column].str.replace(',', '.')
    # If the last column has no name and is all NA, remove it.
    last_column_name = data_df.columns[-1]
    if "Unnamed" in last_column_name and ((data_df.loc[:, last_column_name].isna()) | (data_df.loc[:, last_column_name] == '')).all():
        data_df = data_df.iloc[:, 0:-1]
        
    return data_df




def _find_duplicate_indexes_to_delete(groups, data_df):
    """
    """
    indexes_to_delete = []
    for name, group in groups:
        names = list(data_df.loc[group.index, :].iloc[:, 0])
        fuzz_ratios = compute_fuzz_ratios(names, names)
        fuzz_ratios = {metabolite:series.drop(metabolite) for metabolite, series in fuzz_ratios.items()}
        duplicate_indexes = [group.iloc[i].name for i, (metabolite, series) in enumerate(fuzz_ratios.items()) if len(series) > 0]
        # Add all indexes to list to delete except for 1.
        indexes_to_delete += duplicate_indexes[:-1]
    return indexes_to_delete


def _create_numeric_df(data_df):
    """
    """
    numeric_df = data_df.iloc[:, 1:]
    numeric_df = numeric_df.apply(pandas.to_numeric, errors='coerce')
    # Add metabolites back and drop duplicates.
    numeric_df.loc[:, 'Metabolites'] = data_df.loc[:, 'Metabolite']
    numeric_df = numeric_df.drop_duplicates()
    # Remove Metabolites and find duplicate rows.
    numeric_df = numeric_df.drop('Metabolites', axis=1)
    return numeric_df


