# -*- coding: utf-8 -*-
"""
Metadata Column Matching
~~~~~~~~~~~~~~~~~~~~~~~~
Regular expressions, functions, and classes to match column names and values in mwtab METABOLITES blocks. 

More information can be found on the :doc:`metadata_column_matching` page.
"""

import re

import pandas




# Sometimes 0 is also an NA value, but it can be hard tell unless you see it in an ID column like KEGG ID or something.
# Note the slightly different hyphen character. It is not a duplicate.
NA_VALUES = ['', '-', '−', '--', '---', '.', ',',
             'NA', 'na', 'n.a.', 'N.A.', 'n/a', 'N/A', '<NA>', '#N/A', 'NaN', 'nan', 'N', 'null', 'Null', 'NULL', 'NF',
             'No result', 'NOT Found in Database', 'No ID', 'no data', 'unknown', 'undefined', 'No record', 'NIDB',
             'Not available', 'TBC', 'Internal Standard', 'Intstd', 'internal standard', 'Internal standard',
             'Spiked Stable Isotope Labeled Internal Standards', 'Int Std']

#: Used for wrapping regexes in certain functions.
WRAP_STRING = r'[^a-zA-Z0-9]'
def _create_column_regex_any(match_strings: list[str]) -> str:
    """Return a regular expression string that will match any of the strings in match_strings.
    
    Automatically adds ``WRAP_STRING`` to either side of each string to help with some fuzzy matching. 
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
    
    Automatically adds ``WRAP_STRING`` to either side of each string to help with some fuzzy matching. 
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
    """Used to filter names that match certain criteria.

    Mostly intended to be used through the ColumnFinder class. Created for the purpose 
    of matching tabular column names based on regular expressions and "in" criteria.

    Parameters:
        regex_search_strings: A collection of strings to deliver to re.search() to match a column name. 
          If any string in the collection matches, then the name is matched. This does not simply 
          look for any of the strings within a column name to match. Each string is wrapped with 
          ``WRAP_STRING`` before searching, so the string 'bar' would not be found in 
          the column name "foobarbaz", but would be found in the name "foo bar baz".
        not_regex_search_strings: The same as regex_search_strings, except a match to a column name 
          eliminates that name. Attributes that begin with "not" take precedence over the others. 
          So if a column name matches a string in regex_search_strings and not_regex_search_strings, 
          then it will be filtered OUT.
        regex_search_sets: A collection of sets of strings. Each string in the set of strings must 
          be found in the column name to match, but any set could be found. For example, [('foo', 'bar'), ('baz', 'asd')] 
          will match the name "foo bar" or "bar foo", but not "foobar", due to the aforementioned 
          wrapping with ``WRAP_STRING``. The names "foo", "baz", or "asd" would not match either, but "var asd baz" would.
        in_strings: Similar to regex_search_strings except instead of using re.search() the "in" operator 
          is used. For example, ['foo'] would match the column name "a fool", since 'foo' is in "a fool".
        not_in_strings: The same as in_strings, but matches to a column name eliminate or filter OUT that name.
        in_string_sets: The same as regex_search_sets, but each string in a set is determined to match 
          using the "in" operator instead of re.search(). For example, [('foo', 'bar'), ('baz', 'asd')] 
          WILL match the column name "foobar" because both 'foo' and 'bar' are in the name.
        exact_strings: A collection of strings that must exactly match the column name. For example, 
          ['foo', 'bar'] would only the match the column names "foo" or "bar".

    Examples:
        Find a column for "moverz".
        
        >>> NameMatcher(regex_search_strings = ['m/z', 'mz', 'moverz', 'mx'],
        ...             not_regex_search_strings = ['id'],
        ...             in_strings = ['m.z', 'calcmz', 'medmz', 'm_z', 'obsmz', 'mass to charge', 'mass over z'],
        ...             not_in_strings = ['spec', 'pectrum', 'structure', 'regno', 'retention'])
        
        This is a real example based on the datasets in the Metabolomics Workbench. We can examine 
        some of the strings to illustrate the attributes' function. The "id" string needs to be in 
        not_regex_search_strings, rather than not_in_strings, because "id" is a very small substring 
        that could easily be in a longer word. Putting in not_regex_search_strings means it will 
        most likely match "ID" fields, such as "PubChem ID". Note that it is recommend to lower 
        all column names before filtering and thus use lower case strings, but in general NameMatcher 
        is case sensitive. The "spec" string is in not_in_strings, rather than not_regex_search_strings 
        because the risk of it being in a name that should not be filtered out is low. Also it 
        catches both the full word "spectrum" and its common abbreviation "spec". Hopefully, 
        these 2 explanations of "id" and "spec" have illustrated some of the tradeoffs and 
        advantages of the "in" style attributes versus the "search" style ones.
        
        Find a column for "retention time".
        
        >>> NameMatcher(regex_search_strings = ['rt'],
        ...             regex_search_sets = [['ret', 'time']],
        ...             in_strings = ['rtimes', 'r.t.', 'medrt', 'rtsec', 'bestrt', 'compoundrt', 'rtmed'],
        ...             in_string_sets = [['retention', 'time'], ['rentetion', 'time'], ['retension', 'time']],
        ...             not_in_strings = ['type', 'error', 'index', 'delta', 'feature', 'm/z'])
        
        This is another real example based on the datasets in the Metabolomics Workbench. It 
        illustrates the "set" style attributes quite well. For multi-word column names the 
        "set" style attributes are usually what you want to use. It is possible to to give 
        a string like "retention time", note the space character, to an attribute like "in_strings", 
        but this is more fragile than it seems and won't match some common alternate spellings or 
        mistakes, such as "retention_time" or "retention  time". Using "set" style attributes 
        means you don't have to add as many strings to an attribute like "in_strings". You 
        can still see some repetition in the in_string_sets attribute here though to cover the 
        many mispellings of "retention". "set" style attributes would not be a good use case 
        if the strings in the set must be in a certain order though. The set ['ret', 'time'] 
        will match 'ret' and 'time' in any order. Generally, this will not be a problem because 
        there aren't many instances where you will get a false positive match for a multi-word 
        column due to the order of the words.
        
        Find a column for "other_id".
        
        >>> NameMatcher(not_regex_search_strings = ['cas'],
        ...             in_strings = ['other'],
        ...             in_string_sets = [['database', 'identifier'], ['chemical', 'id'], ['cmpd', 'id'], 
        ...                               ['database', 'id'], ['database', 'match'], ['local', 'id'], 
        ...                               ['row', 'id'], ['comp', 'id'], ['chem', 'id'], ['chro', 'lib', 'id'], 
        ...                               ['lib', 'id']],
        ...             not_in_strings = ['type', 'pubchem', 'chemspider', 'kegg'],
        ...             exact_strings = ['id'],)
        
        This is another real example based on the datasets in the Metabolomics Workbench. It is shown 
        to demonstrate the "exact_strings" attribute. There are many columns that contain the "id" 
        string. There are specific database ID columns, such as those from PubChem or KEGG, but there 
        are often lesser known or individual lab IDs. This example is trying to lump many of the lesser 
        ones into a single "other_id" column. Trying to have "id" in an in_strings or regex_search_strings 
        attribute would cause far too many false positive matches for reasons described in the first 
        example, but there are columns simply labeled "ID", so the only recourse is to use the exact_strings 
        attribute to match them exactly.
        
        Typical usage.
        
        >>> df = pandas.read_csv('some_file.csv')
        >>> name_matcher = NameMatcher(exact_strings = ['foo'])
        >>> modified_columns = {{column_name: column_name.lower().strip() for column_name in df.columns}}
        >>> matching_columns = name_matcher.dict_match(modified_columns)
        
        NameMatcher is really meant to be used as part of a ColumnFinder, but this example uses it 
        directly for simplicity. The instantiated NameMatcher is also very simple in this example 
        because it is trying to show the usage of the dict_match method more than anything else. 
        dict_match requires a dictionary as input, rather than a simple list so that column names 
        can be modified if necessary for easier matching, but then still be linked back to the original 
        name in the dataframe. 

    Attributes:
        regex_search_strings (list[str]): The current list of strings used for regex searching.
        not_regex_search_strings (list[str]): The current list of strings used for regex searching to exclude names.
        regex_search_sets (list[list[str]]): The current list of string sets used for regex searching.
        in_strings (list[str]): The current list of strings used for "in" operator matching.
        not_in_strings (list[str]): The current list of strings used for "in" operator matching to exclude names.
        in_string_sets (list[list[str]]): The current list of string sets used for "in" operator matching.
        exact_strings (list[str]): The current collection of strings used for "==" operator matching.
    """
    
    def __init__(self, regex_search_strings: None|list[str] = None, 
                 not_regex_search_strings: None|list[str] = None, regex_search_sets: None|list[list[str]] = None,
                 in_strings: None|list[str] = None, not_in_strings: None|list[str] = None, 
                 in_string_sets: None|list[list[str]] = None, exact_strings: None|list[str] = None):
        self.regex_search_strings: None|list[str] = regex_search_strings if regex_search_strings else []
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
    """Used to find a mask for certain values in a column.
    
    Mostly intended to be used through the ColumnFinder class. Created for the purpose 
    of matching tabular column data based on regular expressions and type criteria.
    
    Parameters:
        values_type: A string whose only relevant values are 'integer', 'numeric', and 'non-numeric'.
          'integer' will only match values in a column that are integer numbers. 'numeric' will only 
          match values that are numbers, this includes integers. 'non-numeric' will only match values 
          that are non-numeric. Numeric values can be in the value, but cannot be the whole value. 
          For example, '123 id' is considered non-numeric.
        values_regex: A regular expression to positively identify values in a column.
        inverse_values_regex: A regular expression to negatively identify values in a column. This 
          is mutually exclusive with values_regex. If both are given, values_regex takes precedence 
          and inverse_values_regex is ignored. values_type can be combined with either regex and 
          values must match both criteria to match overall.
    
    Examples:
        Simple type example.
        
        >>> vm = ValueMatcher(values_type = 'numeric')
        >>> test = pandas.Series([1, '1', 'foo'])
        >>> vm.series_match(test)
        0     True
        1     True
        2    False
        dtype: bool
        
        This ValueMatcher is very simple and will only match numeric values. Note that numeric 
        values in string form are also recognized as numeric. This is intentional.
        
        Simple regex example.
        
        >>> vm = ValueMatcher(values_regex = 'foo.*')
        >>> test = pandas.Series(['foo', 'bar', 'foobar', 1])
        >>> vm.series_match(test)
        0     True
        1    False
        2     True
        3    False
        dtype: bool
        
        Simple inverse regex example.
        >>> vm = ValueMatcher(values_inverse_regex = 'foo.*')
        >>> test = pandas.Series(['foo', 'bar', 'foobar', 1])
        >>> vm.series_match(test)
        0    False
        1     True
        2    False
        3    False
        dtype: bool
        
        Note that 1 is False in both examples. In general this was designed with strings in mind, so it 
        is recommended to convert all values to strings in any Series delivered to series_match. 
        It is also HIGHLY recommended to use the 'string[pyarrow]' dtype when using the regex attributes. 
        This dtype uses much faster regular expression algorithms and can make orders of magnitude speed 
        differences over Python's built-in regular expressions. There are some features of regular 
        expressions that cannot be used with the 'string[pyarrow]' dtype though. For example, lookahead 
        assertions. More information can be found at https://pypi.org/project/re2/. 
    
    Attributes:
        values_type: The current type of the values being matched.
        values_regex: The regular expression to positively identify values in a column.
        inverse_values_regex: The regular expression to positively exclude values in a column.
    """
    
    def __init__(self, values_type: None|str = None, values_regex: None|str = None, values_inverse_regex: None|str = None):
        self.values_type = values_type if isinstance(values_type, str) else ''
        self.values_regex = values_regex if isinstance(values_regex, str) else ''
        self.values_inverse_regex = values_inverse_regex if isinstance(values_inverse_regex, str) else ''
    
    def series_match(self, series: pandas.Series, na_values: list|None = None, match_na_values: bool = True) -> pandas.Series:
        """Return a mask for the series based on type and regex matching.
        
        "values_regex" and "values_inverse_regex" are mutually exclusive and "values_regex" will take precedence if both are given. 
        "values_type" and one of the regex parameters can both be used, the intermediate masks are ANDed together. 
        "values_type" can only be "integer", "numeric", or "non-numeric" to match those types, respectively.
        
        Args:
            series: series to match values based on type and/or regex.
            na_values: list of values to consider NA values.
            match_na_values: if True, NA values will be consider a match and return True, False otherwise.
        
        Returns:
            A pandas Series the same length as "series" with Boolean values that can be used to select the matching values in the series.
        """
        if na_values is None:
            na_values = []
        
        stripped_series = series.str.strip()
        stripped_series = stripped_series.str.strip('\u200e')
        old_NAs = (stripped_series.isna()) | (stripped_series.isin(na_values))
        
        if self.values_regex:
            regex_match = stripped_series.str.fullmatch(self.values_regex, na=False)
        elif self.values_inverse_regex:
            regex_match = ~stripped_series.str.fullmatch(self.values_inverse_regex, na=True)
        else:
            regex_match = pandas.Series([True]*len(series), index=series.index)
        
        if match_na_values:
            regex_match = regex_match | old_NAs
        
        
        column_to_numeric = pandas.to_numeric(stripped_series, errors='coerce')
        column_to_numeric_NAs = column_to_numeric.isna()
        new_NAs = column_to_numeric_NAs ^ old_NAs
        
        if self.values_type == "integer":
            # The top line will return True for values like '1.0', but the bottom line won't.
            # type_match = (column_to_numeric % 1 == 0) & ~new_NAs
            type_match = ~stripped_series.str.contains('.', regex=False, na=False) & ~new_NAs
        elif self.values_type == "numeric":
            type_match = ~new_NAs
        elif self.values_type == "non-numeric":
            type_match = new_NAs | old_NAs
        else:
            type_match = pandas.Series([True]*len(series), index=series.index)
        
        return regex_match & type_match

class ColumnFinder:
    """Used to find columns in a DataFrame that match a NameMatcher and values in the column that match a ValueMatcher.
    
    This is pretty much just a convenient way to keep the standard_name, NameMatcher, and ValueMatcher 
    together in a single object. Convenience methods to utilize the NameMatcher and ValueMatcher are 
    provided as name_dict_match and values_series_match, respectively.
    
    Parameters:
        standard_name: A string to give a standard name to the column you are trying to find. 
          Not used by any methods.
        name_matcher: The NameMatcher object used to match column names.
        value_matcher: The ValueMatcher object used to match column values.
        
    Examples:
        Basic usage.
        
        >>> df = pandas.DataFrame({'foo':[1, 2, 'asdf'], 'bar':[1, 2, 3]})
        >>> df
            foo  bar
        0     1    1
        1     2    2
        2  asdf    3
        >>> column_finder = ColumnFinder('FOO', NameMatcher(exact_strings = ['foo']), ValueMatcher(values_type = 'numeric'))
        >>> modified_columns = {column_name: column_name.lower().strip() for column_name in df.columns}
        >>> matching_columns = column_finder.name_dict_match(modified_columns)
        >>> matched_column_name = matching_columns[0]
        >>> matched_column_name
        foo
        >>> matching_values = column_finder.values_series_match(df.loc[:, matched_column_name])
        >>> matching_values
        0     True
        1     True
        2    False
        dtype: bool
    
    Attributes:
        standard_name: The standard name of the column trying to be found.
        name_matcher: The NameMatcher object used to match column names.
        value_matcher: The ValueMatcher object used to match column values.
    """
    def __init__(self, standard_name: str, name_matcher: NameMatcher, value_matcher: ValueMatcher):
        self.standard_name = standard_name
        self.name_matcher = name_matcher
        self.value_matcher = value_matcher
    
    def name_dict_match(self, name_map):
        """Convenience method to use the dict_match method for name_matcher.
        """
        return self.name_matcher.dict_match(name_map)
    
    def values_series_match(self, series, na_values = None, match_na_values = True):
        """Convenience method to use the series_match method for value_matcher.
        """
        return self.value_matcher.series_match(series, na_values, match_na_values)
 





def make_list_regex(element_regex: str, delimiter: str , quoted_elements: bool = False, empty_string: bool = False) -> str:
    r"""Creates a regular expression that will match a list of element_regex delimited by delimiter.
    
    Note that delimiter can be a regular expression like (,|;) to match 2 different types of delimiters. 
    If quoted_elements is True, then allow element_regex to be surrounded by single or double quotes. 
    Note that this allows mixed elements, so quoted and unquoted elements are both allowed in the same list. 
    If empty_string is True, then the list regex will match a single element_regex and the empty string. 
    empty_string = True will actually match anything, but the length of the match for strings that 
    are not appropriate will be 0. So this parameter could be useful in some edge case scenarios, 
    but you must investigate the specific match more closely. If the match is the empty string, 
    but the given string is not itself the empty string, then it is not really a match.
    
    Args:
        element_regex: A regular expression in the form of a string that matches the elements of the list to match.
        delimiter: The character(s) that seperate list elements.
        quoted_elements: If True, list elements can be surrounded by single or double quotes.
        empty_string: If True, then allow the returned regular exression to match the empty string.
    
    Returns:
        A regular expression in str form that will match a list of element_regexes delimited by delimiter.
    
    Examples:
        Regular expression to match a list of 4 digit numbers.
        
        >>> regex = make_list_regex(r'\d\d\d\d', r',')
        '((\\d\\d\\d\\d\\s*,\\s*)+(\\d\\d\\d\\d\\s*|\\s*))'
        >>> bool(re.match(regex, '1234'))
        False
        >>> bool(re.match(regex, '1234, 5678'))
        True
        >>> bool(re.match(regex, ''))
        False
               
        Allow the empty string.
        
        >>> regex = make_list_regex(r'\d\d\d\d', r',', empty_string = True)
        '((\\d\\d\\d\\d\\s*,\\s*)*(\\d\\d\\d\\d\\s*|\\s*))'
        >>> bool(re.match(regex, '1234'))
        True
        >>> bool(re.match(regex, '1234, 5678'))
        True
        >>> bool(re.match(regex, ''))
        True
        >>> bool(re.match(regex, 'asdf'))
        True
        >>> re.match(regex, 'asdf')
        <re.Match object; span=(0, 0), match=''>
        
        Allow numbers to be surrounded with quotation marks.
        
        >>> regex = make_list_regex(r'\d\d\d\d', r',', quoted_elements = True)
        '(((\\d\\d\\d\\d\\s*,\\s*)+(\\d\\d\\d\\d\\s*|\\s*))|((\'\\d\\d\\d\\d\'\\s*,\\s*)+(\'\\d\\d\\d\\d\'\\s*|\\s*))|(("\\d\\d\\d\\d"\\s*,\\s*)+("\\d\\d\\d\\d"\\s*|\\s*)))'
        >>> bool(re.match(regex, '1234'))
        False
        >>> bool(re.match(regex, '1234, 5678'))
        True
        >>> bool(re.match(regex, '1234, "5678"'))
        True
        >>> bool(re.match(regex, '\'1234\', "5678"'))
        True
        >>> bool(re.match(regex, ''))
        False
    """
    if quoted_elements:
        # The simplest regex to add quotes around the element, but allows elements like 'element" through.
        # element_regex = r'(\'|"|)' + element_regex + r'(\'|"|)'
        # Forces an element to be surrounded by the same thing on both sides, but still allows for mixed elements.
        # element_regex = r"('" + element_regex + r"'" + '|' + r'"' + element_regex + r'"' + '|' + element_regex + ')'
        
        # All elements must have the same type of quotation marks.
        return '(' + make_list_regex(element_regex, delimiter, empty_string=empty_string) + '|' + \
               make_list_regex(f"'{element_regex}'", delimiter, empty_string=empty_string) + '|' + \
               make_list_regex(f'"{element_regex}"', delimiter, empty_string=empty_string) + ')'
    
    repetition_symbol = '*' if empty_string else '+'
    
    return r'((' + element_regex + r'\s*' + delimiter + r'\s*)' + repetition_symbol + '(' + element_regex + r'\s*|\s*))'



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


# Comment for Sphinx to pull in regular expressions.
INTEGER = r'-?\d+'
FLOAT = r'-?\d*\.\d+'
SCIENTIFIC_NOTATION = r'-?\d*\.\d+E(-|\+)?\d+'
NUMS = '(' + FLOAT + '|' + SCIENTIFIC_NOTATION + '|' + INTEGER + ')'
NUM_RANGE = NUMS + r'(-|\sto\s|−)' + NUMS
LIST_OF_NUMS = make_list_regex(NUMS, ',')
BRACKETED_LIST_OF_NUMS = r'\[' + LIST_OF_NUMS + r'\]'
PARENTHESIZED_LIST_OF_NUMS = r'\(' + LIST_OF_NUMS + r'\)'
LIST_OF_NUMS_UNDERSCORE = make_list_regex(NUMS, '_')
LIST_OF_NUMS_SLASH = make_list_regex(NUMS, '/')
POSITIVE_NUMS = NUMS.replace('-?', '')
POSITIVE_INTS = r'\d+'
LIST_OF_POS_INTS = make_list_regex(POSITIVE_INTS, ',')
LIST_OF_POS_INTS_OR = make_list_regex(POSITIVE_INTS, 'or')
LIST_OF_POS_INTS_BAR = make_list_regex(POSITIVE_INTS, r'\|')
LIST_OF_POS_INTS_SLASH = make_list_regex(POSITIVE_INTS, '/')
LIST_OF_POS_INTS_SEMICOLON = make_list_regex(POSITIVE_INTS, ';')
LIST_OF_POS_INTS_SPACE = make_list_regex(POSITIVE_INTS, ' ')

POSITIVE_FLOATS = r'\d*.\d+'
POSITIVE_SCIENTIFIC_NOTATION = r'\d*\.\d*E(-|\+)?\d+'
POSITIVE_FLOAT_RANGE = POSITIVE_FLOATS + r'\s*(_|-)\s*' + POSITIVE_FLOATS
LIST_OF_POS_FLOATS_UNDERSCORE = make_list_regex(POSITIVE_FLOATS, '_')
POS_FLOAT_PAIRS = '(' +  POSITIVE_FLOATS + r'(_|@)' + POSITIVE_FLOATS + ')'
LIST_OF_POS_FLOAT_PAIRS_UNDERSCORE = make_list_regex(POS_FLOAT_PAIRS, '_')
LIST_OF_POS_FLOAT_PAIRS_NO_SPACE = make_list_regex(POS_FLOAT_PAIRS, '')
LIST_OF_POS_FLOAT_PAIRS_MIXED = make_list_regex(POS_FLOAT_PAIRS, '(//|,)')
POS_INT_FLOAT_PAIR = '(' +  POSITIVE_INTS + r'_' + POSITIVE_FLOATS + ')'

ELEMENT_SYMBOL = (r'([BCFHIKNOPSUVWY]|[ISZ][nr]|[ACELP][ru]|A[cglmst]|B[aehikr]|'
                  r'C[adeflos]|D[bsy]|Es|F[elmr]|G[ade]|H[efgos]|Kr|L[aiv]|M[cdgnot]|'
                  r'N[abdehiop]|O[gs]|P[abdmot]|R[abe-hnu]|S[bcegim]|T[abcehilms]|Xe|Yb)')
ELEMENT_COUNT = r'([1-9]\d*)*'
FORMULA_ELEMENT = ELEMENT_SYMBOL + ELEMENT_COUNT
FORMULA = '(' + FORMULA_ELEMENT + ')+'
LIST_OF_FORMULAS = make_list_regex(FORMULA, ',', True)
BRACKETED_LIST_OF_FORMULAS = r'\[' + LIST_OF_FORMULAS + r'\]'
# C12H22O11, C12 H22 O11, [13C]4H7NO4, [13]C6 H14 [15]N4 O2, [C13]C4H6O5, C21C(13)2H38N7O17P3S, 12C12+14N4+16O19+1H32
ISOTOPIC_NUM = r'\d+'
# Add dueterium.
ISOTOPIC_SYMBOL = ELEMENT_SYMBOL[0:-1] + r'|D)'
ISOTOPIC_ELEMENT = ISOTOPIC_SYMBOL + ELEMENT_COUNT
ISOTOPIC_FORMULA = '(' + ISOTOPIC_ELEMENT + '|' + \
                   r'\[' + ISOTOPIC_NUM + ISOTOPIC_SYMBOL + r'\]' + ELEMENT_COUNT + '|' + \
                   r'\[' + ISOTOPIC_NUM + r'\]' + ISOTOPIC_ELEMENT + '|' + \
                   r'\[' + ISOTOPIC_SYMBOL + ISOTOPIC_NUM + r'\]' + ELEMENT_COUNT + '|' + \
                   ISOTOPIC_SYMBOL + r'\(' + ISOTOPIC_NUM + r'\)' + ELEMENT_COUNT + '|' + \
                   make_list_regex(ISOTOPIC_NUM + ISOTOPIC_SYMBOL + ELEMENT_COUNT, r'\+') + ')+'
LIST_OF_ISOTOPIC_FORMULAS = make_list_regex(ISOTOPIC_FORMULA, ',', True)
BRACKETED_LIST_OF_ISOTOPIC_FORMULAS = r'\[' + LIST_OF_ISOTOPIC_FORMULAS + r'\]'
# C12H22O11+, [C12H22O11]+
CHARGE_FORMULA = r'(\[' + FORMULA + r'\](-|\+)' + '|' + FORMULA + r'(-|\+)' + ')'
# CH3(CH2)16COOH
GROUP_FORMULA = '(' + FORMULA_ELEMENT + '|' + r'\(' + FORMULA + r'\)\d+' + r')+'

ORGANIC_ELEMENT_SYMBOL = r'([CHNOPS])'
ORGANIC_FORMULA_ELEMENT = ORGANIC_ELEMENT_SYMBOL + ELEMENT_COUNT
ORGANIC_FORMULA = '(' + ORGANIC_FORMULA_ELEMENT + '){4,}'

SMILES_ELEMENT_SYMBOL = r'\d?[a-zA-Z][a-z]?'
# Note the , , , and  symbols at the end of the character set are '\x01', 'x02', '\x03', and '\x04'. 
# I don't think they are apart of the SMILES characters, but they appeared in the SMILES of some datasets.
# There are some SMILES in AN003143 like 'C[n]1c[n]cc1C[C@@H](NC(=O)CCN)C(O)=O |&1:7|'. 
# I don't think the end bit between '|' is part of a legit SMILES, but I added it to pass this dataset.
SMILES = r'(\[?' + SMILES_ELEMENT_SYMBOL + r'[0-9@+\-[\]()\\/%=#$.*]*)+' + r'( \|[0-9:&,w]+\|)?'
LIST_OF_SMILES_SEMICOLON = make_list_regex(SMILES, ';')

INCHIKEY = r'(InChIKey=)?[a-zA-Z]{14}-[a-zA-Z]{10}-[a-zA-Z]?'
INCHIKEY_OR_NULL = '(' + INCHIKEY + r'|null|No record)'
LIST_OF_INCHIKEYS = '(' + INCHIKEY + r'\s*,\s*)+' + '(' + INCHIKEY + r'\s*|\s*)'
LIST_OF_INCHIKEYS_SLASH = LIST_OF_INCHIKEYS.replace(',', '/')

# InChI=1 or InChI=1 where "1" is a version number and "S" means it's a standard InChI.
INCHI_HEAD = r'InChI='
INCHI_VERSION = r'\d+S?'
INCHI_FORMULA = r'/' + FORMULA
INCHI_SKELETAL_LAYER = r'/c(\d+([-,()])?)+'
INCHI_HYDROGEN_LAYER = r'/h' + '(' + make_list_regex(r'\d+(-\d+)?H?\d*', ',') + r')?' + r'(\(H\d*-?(,\d+)+\))*'
INCHI_CHARGE_LAYER = r'/q(-|\+)\d+'
INCHI_PROTONATION_LAYER = r'/p(-|\+)\d+'
INCHI_STEREOCHEMISTRY_LAYER = r'/t' + make_list_regex(r'(\d+(-|\+|\?|u)|M)', ',')
INCHI_STEREOCHEMISTRY_SUBLAYER1 = r'/m\d+'
INCHI_STEREOCHEMISTRY_SUBLAYER2 = r'/s\d+'
INCHI_FIXED_HYDROGEN_LAYER = r'/f(' + FORMULA + r')?' + \
                             '(' + INCHI_HYDROGEN_LAYER + r')?' + \
                             '(' + INCHI_STEREOCHEMISTRY_SUBLAYER2 + r')?' + \
                             '(' + INCHI_CHARGE_LAYER + r')?'
INCHI_DOUBLEBOND_LAYER = r'/b' + make_list_regex(r'\d+(-|\+)\d+(-|\+)', ',')
INCHI_ISOTOPIC_LAYER = r'/i' + make_list_regex(r'(\d+(-|\+)\d+|\d+[A-Z]\d*|\d+(-|\+)\d+[A-Z]?\d*)', ',') + '(' + INCHI_STEREOCHEMISTRY_SUBLAYER2 + r')?'

FULL_INCHI = '(' + INCHI_HEAD + r')?' + \
             INCHI_VERSION + \
             '(' + INCHI_FORMULA + r')?' + \
             '(' + INCHI_SKELETAL_LAYER + r')?' + \
             '(' + INCHI_HYDROGEN_LAYER + r')?' + \
             '(' + INCHI_CHARGE_LAYER + r')?' + \
             '(' + INCHI_PROTONATION_LAYER + r')?' + \
             '(' + INCHI_DOUBLEBOND_LAYER + r')?' + \
             '(' + INCHI_STEREOCHEMISTRY_LAYER + r')?' + \
             '(' + INCHI_STEREOCHEMISTRY_SUBLAYER1 + r')?' + \
             '(' + INCHI_STEREOCHEMISTRY_SUBLAYER2 + r')?' + \
             '(' + INCHI_ISOTOPIC_LAYER + r')?' + \
             '(' + INCHI_FIXED_HYDROGEN_LAYER + r')?'

CHEAP_INCHI = r'\s*' + '(' + INCHI_HEAD + r')?' + INCHI_VERSION + r'((/c|/h|/i|/t|/m|/s|/q|/f|/p|/b|/[A-Z])(\S)+)+' + r'\s*'
LIST_OF_INCHI = make_list_regex(CHEAP_INCHI, ',', True)
BRACKETED_LIST_OF_INCHI = r'\[' + LIST_OF_INCHI + r'\]'

# M+H, [M+H]+, -H(-), M+AGN+H, M+Acid, (M-H)/2, [M+NH4] +_[M+Na]+, [M+H–H2O]+, [M-2H](2-)
# Cat, Hac, and Chol-head are wierd special cases.
ION_ELEMENTS = r'(\d?(m|M)' + '|' + \
               r'(-|\+)?\d*' + FORMULA + r'\d*(\((-|\+)\))?' +  '|' + \
               r'[ [a-zA-Z]*(A|a)cid' + '|' + \
               r'Cat' + '|' + r'Chol-head' + '|' + r'Hac' + '|' + r'\di' + '|' + r'FA|NA|A' + ')'
ION = '(' + make_list_regex(ION_ELEMENTS, r'(-|\+)?') + ')'
ION = '(' + ION + '|' + \
      r'\[' + ION + r'\]\(?\d?\s?(-|\+)?\d?\)?' + '|' + \
      ION + r'\](-|\+)?' + '|' + \
      r'\(' + ION + r'\)((/\d)|(-|\+))?' + ')'
LIST_OF_IONS = make_list_regex(ION, ',')
LIST_OF_IONS_SPACE = make_list_regex(ION, ' ')
LIST_OF_IONS_UNDERSCORE = make_list_regex(ION, '_')
LIST_OF_IONS_MIXED = make_list_regex(ION, '(_| )')
LIST_OF_IONS_NO_DELIMITER = make_list_regex(ION, '')
BRACKETED_LIST_OF_IONS = r'\[' + make_list_regex(ION, ',', True) + r'\]'

# There are values like CA1511 that I can't confirm are KEGG values. You can't find them in the compound database, but they appear often.
# The prefixes CE, UP, Z, and U are the same. There is a common mistake of only having 4 numbers after C, so that is accounted for as well.
KEGG = '(' + r'(cpd:)?[CDMGKRZU]0?\d{5}\?{0,2}' + '|' + r'(DG|ko)\d{5}' + '|' + r'(CA|CE|UP|C)\d{4}' + '|' + r'(NA|n/a)' + ')'
LIST_OF_KEGG = make_list_regex(KEGG, ',')
LIST_OF_KEGG_SEMICOLON = make_list_regex(KEGG, ';')
LIST_OF_KEGG_SLASH = make_list_regex(KEGG, '/')
LIST_OF_KEGG_SPACE = make_list_regex(KEGG, ' ')
LIST_OF_KEGG_DOUBLE_SLASH = make_list_regex(KEGG, '//')
LIST_OF_KEGG_UNDERSCORE = make_list_regex(KEGG, '_')
LIST_OF_KEGG_HYPHEN = make_list_regex(KEGG, '-')
LIST_OF_KEGG_MIXED = make_list_regex(KEGG, '(/|,)')
LIST_OF_KEGG_BAR = make_list_regex(KEGG, r'(\|)')

HMDB = '(' + r'(HMDB|HDMB|YMDB|HMBD)\d+(\*|\?)?' + '|' + r'n/a' + ')'
LIST_OF_HMDB = make_list_regex(HMDB, ',')
LIST_OF_HMDB_SLASH = make_list_regex(HMDB, '/')
LIST_OF_HMDB_AMPERSAND = make_list_regex(HMDB, '&')
LIST_OF_HMDB_SEMICOLON = make_list_regex(HMDB, ';')
LIST_OF_HMDB_SPACE = make_list_regex(HMDB, ' ')
LIST_OF_HMDB_UNDERSCORE = make_list_regex(HMDB, '_')
HMDB_INT = r'\d{,5}'
LIST_OF_HMDB_INTS = make_list_regex(HMDB_INT, ',')
LIST_OF_HMDB_INTS_SLASH = make_list_regex(HMDB_INT, '/')

LIPID_MAPS = '(' + r'LM(PK|ST|GL|FA|SP|GP|PR|SL)[0-9A-Z]{8,10}\*?' + '|' + r'(ST|FA|PR|GP|PK|GL|SP)\d{4,6}-' + FORMULA + ')'
LIST_OF_LMP = make_list_regex(LIPID_MAPS, ',')
LIST_OF_LMP_UNDERSCORE = make_list_regex(LIPID_MAPS, '_')
LIST_OF_LMP_SLASH = make_list_regex(LIPID_MAPS, '/')

# When using Excel, a CAS number can get mistaken for a date and it will automatically change the value.
DATE = r'\d{1,2}/\d{1,2}/(\d{4}|\d{2})'
CAS = r'(CAS: ?)?\d+-\d\d-0?\d' + '|' + DATE
LIST_OF_CAS = make_list_regex(CAS, ',')
LIST_OF_CAS_SEMICOLON = make_list_regex(CAS, ';')
# Comment for Sphinx to find the end of regular expressions.


# Important NOTE! All of the regular expressions delivered to ValueMatchers are surrounded with an additional set of 
# parenthesis. This is because they are later used with the pyarrow integration in pandas, and the behavior is 
# slightly different between python's built in regular expressions and pyarrow's. TLDR, add wrapping parenthesis 
# to avoid problems matching using pyarrow.
column_finders = [
    ColumnFinder("moverz_quant",
                 NameMatcher(regex_search_strings = ['m/z', 'mz', 'moverz', 'mx'],
                             not_regex_search_strings = ['id'],
                             in_strings = ['m.z', 'calcmz', 'medmz', 'm_z', 'obsmz', 'mass to charge', 'mass over z'],
                             not_in_strings = ['spec', 'pectrum', 'structure', 'regno', 'retention'],),
                 ValueMatcher(values_regex = '(' + NUMS + '|' + \
                                             LIST_OF_NUMS + '|' + \
                                             LIST_OF_NUMS_UNDERSCORE + '|' + \
                                             NUMS + r'\s*/\s*' + NUMS + '|' + \
                                             '(' + NUMS + r'\s*\(' + NUMS + r'\)' + r'\s*;\s*)+' + '(' + NUMS + r'\s*\(' + NUMS + r'\)' + r'\s*|\s*)' + '|' + \
                                             NUMS + r'(\s*>\s*|\s*<\s*)' + NUMS + '|' + \
                                             '(' + NUMS + r'\s*)+' + '|' + \
                                             NUMS + r'\s*\(\s*' + NUMS + r'\s*\)' + '|' + \
                                             NUMS + r'\s*-\s*' + NUMS + ')',)),
    
    ColumnFinder("mass",
                 NameMatcher(regex_search_strings = ['mass', 'quantmass', 'masses', 'mw', 'weight'],
                             not_regex_search_strings = ['id'],
                             in_strings = ['exactmass', 'obsmass', 'calcmass', 'monoisotopicmass', 'molwt'],
                             not_in_strings = ['spec', 'pectrum', 'structure', 'regno', 'charge', 'over z', 'rsd', 'm/z'],
                             exact_strings = ['m meas.'],),
                 ValueMatcher(values_regex = '(' + NUMS + '|' + \
                                             LIST_OF_NUMS + '|' + \
                                             NUMS + r'\s*/\s*' + NUMS + r'(\s*Da)?|' + \
                                             NUMS + r'\s*-\s*' + NUMS + ')',)),
    
    ColumnFinder("parent_moverz_quant",
                 NameMatcher(exact_strings = ['parent'],),
                 ValueMatcher(values_type = 'numeric',)),
    
    ColumnFinder("mass_spectrum",
                 NameMatcher(in_strings = ['spec', 'pectrum'],
                             not_in_strings = ['species', 'composite'],),
                 ValueMatcher(values_regex = '(' + '(' + NUMS + ':' + NUMS + r'(_|\s+|$))+' + '|' + \
                                             NUMS + '|' + \
                                             LIST_OF_NUMS + ')',)),
    
    ColumnFinder("composite_mass_spectrum",
                 NameMatcher(in_string_sets = [['composite', 'spectrum']],
                             not_in_strings = ['species'],),
                 ValueMatcher(values_regex = '(' + '(' + PARENTHESIZED_LIST_OF_NUMS + r'\s*)+' + ')',)),
    
    ColumnFinder("inchi_key",
                 NameMatcher(in_strings = ['inchikey', 'inchi-key', 'inchi_key', 'inchi key'],),
                 ValueMatcher(values_regex = '(' + INCHIKEY + r'(\*?|\?*)' + '|' +\
                                             INCHIKEY_OR_NULL + r'(\s*or\s*|\s*;\s*|\s*_\s*|\s*&\s*)' + INCHIKEY_OR_NULL + '|' +\
                                             r'Sum\s*\(\s*' + INCHIKEY + r'(\s*\+\s*)' + INCHIKEY + r'\s*\)' + '|' +\
                                             LIST_OF_INCHIKEYS_SLASH + ')',)),
    
    ColumnFinder("inchi",
                 NameMatcher(in_strings = ['inchi'],
                             not_in_strings = ['key'],),
                 ValueMatcher(values_regex = '(' + CHEAP_INCHI + '|' + BRACKETED_LIST_OF_INCHI + ')',)),
    
    ColumnFinder("smiles",
                 NameMatcher(in_strings = ['smile'],),
                 ValueMatcher(values_regex = '(' + SMILES + '|' + LIST_OF_SMILES_SEMICOLON + ')',)),
    
    ColumnFinder("formula",
                 NameMatcher(in_strings = ['formula'],),
                 ValueMatcher(values_type = 'non-numeric',
                              values_regex = '(' + ISOTOPIC_FORMULA.replace(r'\d*', r'\d*\s*') + '|' + \
                                             CHARGE_FORMULA + '|' + \
                                             GROUP_FORMULA + '|' + \
                                             BRACKETED_LIST_OF_FORMULAS + '|' + \
                                             BRACKETED_LIST_OF_ISOTOPIC_FORMULAS + ')',)),
    
    # Metabolite is the first left-most column in every table for mwtab. This Finder will only create false positives.
    # ColumnFinder("metabolite",
    #              NameMatcher(exact_strings = ['metabolite'],),
    #              ValueMatcher(values_type = 'non-numeric',)),
    
    ColumnFinder("compound",
                 NameMatcher(in_strings = ['compound', 'compund'],
                             not_in_strings = ['kegg', 'formula', 'pubchem', 'mass', 'rt', 'algo', 'id', 'name'],),
                 ValueMatcher(values_type = 'non-numeric',)),
    
    ColumnFinder("name",
                 NameMatcher(in_strings = ['name'],
                             in_string_sets = [['name', 'refmet']],
                             not_in_strings = ['adduct', 'named', 'internal', 'ion', 'metabolite_name'],),
                 ValueMatcher(values_type = 'non-numeric',)),
    
    ColumnFinder("refmet",
                 NameMatcher(in_strings = ['refmet'],
                             not_in_strings = ['name', 'in'],),
                 ValueMatcher(values_regex = '(' + POSITIVE_INTS + ')',)),
    
    ColumnFinder("class",
                 NameMatcher(in_strings = ['class']),
                 ValueMatcher(values_type = 'non-numeric',
                              values_inverse_regex = '(' + ORGANIC_FORMULA + ')',)),
    
    ColumnFinder("pathway",
                 NameMatcher(in_strings = ['pathway'],
                             not_in_strings = ['sort'],),
                 ValueMatcher(values_type = 'non-numeric',)),
    
    ColumnFinder("pathway_sortorder",
                 NameMatcher(in_string_sets = [['pathway', 'sort']],),
                 ValueMatcher(values_type = 'integer',)),
    
    ColumnFinder("ion",
                 NameMatcher(regex_search_strings = ['ion', 'ions'],
                             not_in_strings = ['adduct', 'm/z', 'mass'],),
                 ValueMatcher(values_regex = '(' + ION + '|' + \
                                             LIST_OF_IONS + '|' + \
                                             LIST_OF_IONS_SPACE + '|' + \
                                             NUMS + '|' + \
                                             LIST_OF_NUMS + '|' + \
                                             NUMS + r'(\s*>\s*|\s*<\s*)' + NUMS + '|' + \
                                             LIST_OF_NUMS_SLASH + ')',)),
    
    ColumnFinder("adduct",
                 NameMatcher(in_strings = ['adduct'],
                             not_in_strings = ['formula'],),
                 ValueMatcher(values_regex = '(' + ION + '|' + \
                                             LIST_OF_IONS + '|' + \
                                             LIST_OF_IONS_SPACE + '|' + \
                                             LIST_OF_IONS_UNDERSCORE + '|' + \
                                             LIST_OF_IONS_MIXED + '|' + \
                                             LIST_OF_IONS_NO_DELIMITER + '|' + \
                                             BRACKETED_LIST_OF_IONS + ')',)),
    
    ColumnFinder("species",
                 NameMatcher(in_strings = ['species'],
                             not_in_strings = ['is_species', 'ion'],),
                 ValueMatcher(values_regex = '(' + ION + '|' + \
                                             LIST_OF_IONS + '|' + \
                                             LIST_OF_IONS_UNDERSCORE + ')',)),
    
    ColumnFinder("pubchem_id",
                 NameMatcher(regex_search_strings = ['cid'],
                             in_strings = ['pubchem'],
                             not_in_strings = ['formula', 'kegg'],),
                 ValueMatcher(values_regex = '(' + POSITIVE_INTS + r'[&?]?' + '|' + \
                                             LIST_OF_POS_INTS + '|' + \
                                             LIST_OF_POS_INTS_OR + '|' + \
                                             LIST_OF_POS_INTS_SLASH + '|' + \
                                             LIST_OF_POS_INTS_SPACE + '|' + \
                                             LIST_OF_POS_INTS_SEMICOLON + '|' + \
                                             r'Sum \(\d+ \+ \d+\)' + '|' + \
                                             r'CID' + POSITIVE_INTS + ')',)),
    
    ColumnFinder("kegg_id",
                 NameMatcher(in_strings = ['kegg'],
                             not_in_strings = ['name'],),
                 ValueMatcher(values_regex = '(' + KEGG + '|' + \
                                             LIST_OF_KEGG + '|' + \
                                             LIST_OF_KEGG_SEMICOLON + '|' + \
                                             LIST_OF_KEGG_SLASH + '|' + \
                                             LIST_OF_KEGG_DOUBLE_SLASH + '|' + \
                                             LIST_OF_KEGG_UNDERSCORE + '|' + \
                                             LIST_OF_KEGG_HYPHEN + '|' + \
                                             LIST_OF_KEGG_MIXED + '|' + \
                                             LIST_OF_KEGG_SPACE + '|' + \
                                             LIST_OF_KEGG_BAR + '|' + \
                                             KEGG + r'-' + FORMULA  + '|' + \
                                             KEGG + r';\d+' + ')',)),
    
    ColumnFinder("hmdb_id",
                 NameMatcher(in_strings = ['hmdb', 'human metabolome'],
                             in_string_sets = [['hmp', 'id']],
                             not_in_strings = ['class'],),
                 ValueMatcher(values_regex = '(' + HMDB + '|' + \
                                             HMDB_INT + '|' + \
                                             LIST_OF_HMDB + '|' + \
                                             LIST_OF_HMDB_SLASH + '|' + \
                                             LIST_OF_HMDB_INTS + '|' + \
                                             LIST_OF_HMDB_INTS_SLASH + '|' + \
                                             r'Sum \(HMDB\d+ \+ HMDB\d+\)' + '|' + \
                                             LIST_OF_HMDB_AMPERSAND + '|' + \
                                             LIST_OF_HMDB_SEMICOLON + '|' + \
                                             LIST_OF_HMDB_SPACE + '|' + \
                                             r'METPA\d+' + '|' + \
                                             LIST_OF_HMDB_UNDERSCORE + ')',)),
    
    ColumnFinder("lm_id",
                 NameMatcher(in_strings = ['lipidmaps', 'lmid'],
                             in_string_sets = [['lmp', 'id'], ['lipid', 'map'], ['lm', 'id']],),
                 ValueMatcher(values_regex = '(' + LIPID_MAPS + '|' + \
                                             LIST_OF_LMP_UNDERSCORE + '|' + \
                                             LIST_OF_LMP + '|' + \
                                             LIST_OF_LMP_SLASH + '|' + \
                                             POSITIVE_INTS + ')',)),
    
    ColumnFinder("chemspider_id",
                 NameMatcher(in_strings = ['chemspider'],),
                 ValueMatcher(values_regex = '(' + POSITIVE_INTS + r'[&?]?' + '|' + \
                                             LIST_OF_POS_INTS + '|' + \
                                             LIST_OF_POS_INTS_OR + '|' + \
                                             LIST_OF_POS_INTS_SLASH + '|' + \
                                             LIST_OF_POS_INTS_SPACE + '|' + \
                                             LIST_OF_POS_INTS_SEMICOLON + '|' + \
                                             r'Sum \(\d+ \+ \d+\)' + '|' + \
                                             'CID' + POSITIVE_INTS + '|' + \
                                             'CSID' + POSITIVE_INTS + '|' + \
                                             r'[CD]\d{5}' + '|' + \
                                             make_list_regex(r'[CD]\d{5}', r'\|') + ')',)),
    
    ColumnFinder("metlin_id",
                 NameMatcher(in_strings = ['metlin'],),
                 ValueMatcher(values_regex = '(' + POSITIVE_INTS + '|' + r'METLIN:' + POSITIVE_INTS + ')',)),
    
    ColumnFinder("cas_number",
                 NameMatcher(regex_search_strings = ['cas'],),
                 ValueMatcher(values_regex = '(' + CAS + '|' + LIST_OF_CAS + '|' + LIST_OF_CAS_SEMICOLON + ')',)),
    
    ColumnFinder("binbase_id",
                 NameMatcher(regex_search_strings = ['bb'],
                             in_strings = ['binbase'],),
                 ValueMatcher(values_regex = '(' + POSITIVE_INTS + ')',)),
    
    ColumnFinder("chebi_id",
                 NameMatcher(in_strings = ['chebi'],),
                 ValueMatcher(values_regex = '(' + POSITIVE_INTS + ')',)),
    
    ColumnFinder("mw_regno",
                 NameMatcher(in_strings = ['regno'],
                             in_string_sets = [['mw', 'structure']],),
                 ValueMatcher(values_regex = '(' + POSITIVE_INTS + '|' + LIST_OF_POS_INTS + ')',)),
    
    ColumnFinder("mzcloud_id",
                 NameMatcher(in_string_sets = [['mz', 'cloud', 'id']],),
                 ValueMatcher(values_regex = '(' + r'((Reference|Autoprocessed)-)?' + POSITIVE_INTS + ')',)),
    
    ColumnFinder("identifier",
                 NameMatcher(in_strings = ['identifier', 'retention time_m/z', 'feature@rt'],
                             not_in_strings = ['pubchem', 'study', 'database'],),
                 ValueMatcher(values_regex = '(' + POS_FLOAT_PAIRS + r'(n|m/z)?' + '|' + \
                                             POS_INT_FLOAT_PAIR + '|' + \
                                             LIST_OF_POS_FLOAT_PAIRS_UNDERSCORE + '|' + \
                                             LIST_OF_POS_FLOAT_PAIRS_NO_SPACE + '|' + \
                                             LIST_OF_POS_FLOAT_PAIRS_MIXED + '|' + \
                                             r'CHEBI:\d+' + ')',)),
    
    ColumnFinder("other_id",
                 NameMatcher(not_regex_search_strings = ['cas'],
                             in_strings = ['other'],
                             in_string_sets = [['database', 'identifier'], ['chemical', 'id'], ['cmpd', 'id'], ['database', 'id'], ['database', 'match'], ['local', 'id'], ['row', 'id'], ['comp', 'id'], ['chem', 'id'], ['chro', 'lib', 'id'], ['lib', 'id']],
                             not_in_strings = ['type', 'pubchem', 'chemspider', 'kegg'],
                             exact_strings = ['id'],),
                 ValueMatcher(values_inverse_regex = '(' + FLOAT + '|' + SCIENTIFIC_NOTATION + ')',)),
    
    ColumnFinder("other_id_type",
                 NameMatcher(in_string_sets = [['other', 'type'], ['source', 'database']],),
                 ValueMatcher(values_type = 'non-numeric',)),
    
    ColumnFinder("retention_time",
                 NameMatcher(regex_search_strings = ['rt'],
                             regex_search_sets = [['ret', 'time']],
                             in_strings = ['rtimes', 'r.t.', 'medrt', 'rtsec', 'bestrt', 'compoundrt', 'rtmed'],
                             in_string_sets = [['retention', 'time'], ['rentetion', 'time'], ['retension', 'time']],
                             not_in_strings = ['type', 'error', 'index', 'delta', 'feature', 'm/z'],),
                 ValueMatcher(values_regex = '(' + POSITIVE_FLOATS + '|' + \
                                             r'\d' + '|' + \
                                             POSITIVE_SCIENTIFIC_NOTATION + '|' + \
                                             POSITIVE_FLOAT_RANGE + '|' + \
                                             LIST_OF_POS_FLOATS_UNDERSCORE + ')',)),
    
    ColumnFinder("delta_rt",
                 NameMatcher(in_strings = ['deltart'],
                             in_string_sets = [['delta', 'rt']],
                             not_in_strings = ['type', 'error', 'index'],),
                 ValueMatcher(values_regex = FLOAT + '|' + r'0',)),
    
    ColumnFinder("retention_index",
                 NameMatcher(regex_search_strings = ['ri'],
                             regex_search_sets = [['ret', 'ind'], ['ret', 'index']],
                             in_strings = ['rindex'],
                             in_string_sets = [['retention', 'index'], ['rentetion', 'index'], ['reten', 'index']],
                             not_in_strings = ['type', 'error'],),
                 ValueMatcher(values_regex = '(' + POSITIVE_FLOATS + '|' + \
                                             r'\d' + '|' + \
                                             POSITIVE_SCIENTIFIC_NOTATION + '|' + \
                                             POSITIVE_FLOAT_RANGE + '|' + \
                                             LIST_OF_POS_FLOATS_UNDERSCORE + ')',)),
    
    ColumnFinder("retention_index_type",
                 NameMatcher(in_string_sets = [['retention', 'index', 'type'], ['ri', 'type']],
                             not_in_strings = ['error'],),
                 ValueMatcher(values_type = 'non-numeric',)),
    
    ColumnFinder("abbreviation",
                 NameMatcher(in_strings = ['abbreviation'],),
                 ValueMatcher(values_type = 'non-numeric',)),
    
    ColumnFinder("assignment_certainty",
                 NameMatcher(in_string_sets = [['assignment', 'certainty']],),
                 ValueMatcher(values_regex = '(' + POSITIVE_INTS + ')',)),
    
    ColumnFinder("comment",
                 NameMatcher(in_strings = ['comment'],),
                 ValueMatcher(values_inverse_regex = '(' + FLOAT + '|' + SCIENTIFIC_NOTATION + '|' + r'\d{2,}' + ')',)),
    
    ColumnFinder("assignment_method",
                 NameMatcher(in_string_sets = [['assignment', 'method']],),
                 ValueMatcher(values_type = 'non-numeric',)),
    
    ColumnFinder("isotopologue",
                 NameMatcher(in_strings = ['isotopologue'],
                             in_string_sets = [['isotope', 'count']],
                             not_in_strings = ['type'],),
                 ValueMatcher(values_type = 'non-numeric',)),
    
    ColumnFinder("isotopologue_type",
                 NameMatcher(in_strings = ['isotopologue%type', 'isotope'],
                             not_in_strings = ['count'],),
                 ValueMatcher(values_type = 'non-numeric',)),
    
    ColumnFinder("peak_description",
                 NameMatcher(in_string_sets = [['peak', 'description']],),
                 ValueMatcher(values_type = 'non-numeric',)),
    
    ColumnFinder("peak_pattern",
                 NameMatcher(in_string_sets = [['peak', 'pattern']],),
                 ValueMatcher(values_type = 'non-numeric',)),
    
    ColumnFinder("transient_peak",
                 NameMatcher(in_string_sets = [['transient', 'peak']],
                             not_in_strings = ['type'],),
                 ValueMatcher(values_regex = '(' + POSITIVE_INTS + ')',)),
    
    ColumnFinder("transient_peak_type",
                 NameMatcher(in_string_sets = [['transient', 'peak', 'type']],),
                 ValueMatcher(values_type = 'non-numeric',)),
    
    ColumnFinder("fish_coverage",
                 NameMatcher(in_string_sets = [['fish', 'coverage']],),
                 ValueMatcher(values_regex = '(' + POSITIVE_NUMS + ')',)),
    
    ColumnFinder("msi_category",
                 NameMatcher(in_strings = ['msicategory'],),
                 ValueMatcher(values_regex = r'1',)),
    
    ColumnFinder("annotations",
                 NameMatcher(not_regex_search_strings = ['id'],
                             in_strings = ['annotation'],
                             not_in_strings = ['source', 'approach', 'confidence', 'level'],),
                 ValueMatcher(values_type = 'non-numeric',)),
    
    ColumnFinder("istd",
                 NameMatcher(in_strings = ['internal'],
                             exact_strings = ['istd'],),
                 ValueMatcher(values_type = 'non-numeric',)),
    
    ColumnFinder("platform",
                 NameMatcher(in_strings = ['platform'],),
                 ValueMatcher(values_type = 'non-numeric',)),
    
    ColumnFinder("ms_method",
                 NameMatcher(in_strings = ['method'],
                             not_in_strings = ['assignment'],),
                 ValueMatcher(values_type = 'non-numeric',)),
    
    ColumnFinder("polarity",
                 NameMatcher(in_strings = ['polarity'],),
                 ValueMatcher(values_regex = r'(?i)((neg|pos|1|-1|\+|positive|negative)' + r'|\[M\+H\]\+|\[M-H\]-|5MM\+|5MM-)',)),
    
    ColumnFinder("esi_mode",
                 NameMatcher(in_strings = ['esi'],),
                 ValueMatcher(values_regex = r'(neg|pos|1|-1|(ESI )?\(\+\)( ESI)?|(ESI )?\(-\)( ESI| ES\))?|positive|negative)',)),
    
    ColumnFinder("ionization_mode",
                 NameMatcher(regex_search_sets = [['pos', 'neg']],
                             in_strings = ['ionization', 'ionisation'],
                             not_in_strings = ['confirmed'],
                             exact_strings = ['mode', 'ms mode'],),
                 ValueMatcher(values_regex = r'(?i)((neg|pos|1|-1|(ES)?\+|(ES)?-|positive|negative|TOF|Splitless|Split30))',)),
    
    ColumnFinder("frequency",
                 NameMatcher(in_strings = ['frequency'],),
                 ValueMatcher(values_regex = '(' + POSITIVE_INTS + ')',)),
]

column_finders = {finder.standard_name: finder for finder in column_finders}


implied_pairs = {'other_id': ['other_id_type'], 'retention_index': ['retention_index_type']}


