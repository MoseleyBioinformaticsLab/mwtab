#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
mwtab.tokenizer
~~~~~~~~~~~~~~~

This module provides the :func:`~mwtab.tokenizer.tokenizer` lexical analyzer for
`mwTab` format syntax. It is implemented as Python generator-based state
machine which generates (yields) tokens one at a time when :py:func:`next()`
is invoked on :func:`~mwtab.tokenizer.tokenizer` instance.

Each token is a tuple of "key-value"-like pairs, tuple of
``SUBJECT_SAMPLE_FACTORS`` or tuple of data deposited between
``*_START`` and ``*_END`` blocks.
"""

from __future__ import print_function, division, unicode_literals
from collections import deque, namedtuple
import re


KeyValue = namedtuple("KeyValue", ["key", "value"])

def _results_file_line_to_dict(line: str):
    """Parse a RESULTS_FILE line into a dictionary.
    
    Args:
        line: The line to parse. Expected to just be the value without "RESULTS_FILE" in it.
    
    Returns:
        A dictionary of the values found in the line, won't always have every key.
    """
    filename_regex = '(.*?)([^\s]+?\s*?)((\s(UNITS|Has m/z|Has RT|RT units))|$)'
    units_regex = '(.*)UNITS:(.*?)((\s(Has m/z|Has RT|RT units))|$)'
    has_mz_regex = '(.*)Has m/z:(.*?)((\s(UNITS|Has RT|RT units))|$)'
    has_rt_regex = '(.*)Has RT:(.*?)((\s(Has m/z|UNITS|RT units))|$)'
    rt_units_regex = '(.*)RT units:(.*?)((\s(Has m/z|Has RT|UNITS))|$)'
    
    filename_value = None if not (match := re.match(filename_regex, line)) else match.group(2).rstrip()
    units_value = None if not (match := re.match(units_regex, line)) else match.group(2).strip()
    has_mz_value = None if not (match := re.match(has_mz_regex, line)) else match.group(2).strip()
    has_rt_value = None if not (match := re.match(has_rt_regex, line)) else match.group(2).strip()
    rt_units_value = None if not (match := re.match(rt_units_regex, line)) else match.group(2).strip()
    
    results_file_dict = {}
    if filename_value:
        results_file_dict['filename'] = filename_value
    if units_value:
        results_file_dict['UNITS'] = units_value
    if has_mz_value:
        results_file_dict['Has m/z'] = has_mz_value
    if has_rt_value:
        results_file_dict['Has RT'] = has_rt_value
    if rt_units_value:
        results_file_dict['RT units'] = rt_units_value
    
    return results_file_dict


def tokenizer(text, dict_type = None):
    """A lexical analyzer for the `mwtab` formatted files.

    :param text: `mwTab` formatted text.
    :type text: py:class:`str`
    :param dict_type: the type of dictionary to use, default is dict.
    :return: Tuples of data.
    :rtype: py:class:`~collections.namedtuple`
    """
    if dict_type is None:
        dict_type = dict
        
    stream = deque(text.split("\n"))

    while len(stream) > 0:
        line = stream.popleft()
        try:

            # header
            if line.startswith("#METABOLOMICS WORKBENCH"):
                yield KeyValue("#METABOLOMICS WORKBENCH", "\n")
                for i, identifier in enumerate(re.split(r' |\t', line)):
                    if ":" in identifier:
                        key, value = identifier.split(":")
                        yield KeyValue(key, value)

            # SUBJECT_SAMPLE_FACTORS header (reached new section)
            elif line.startswith("#SUBJECT_SAMPLE_FACTORS:"):
                yield KeyValue("#ENDSECTION", "\n")
                yield KeyValue("#SUBJECT_SAMPLE_FACTORS", "\n")

            # section header (reached new section)
            elif line.startswith("#"):
                yield KeyValue("#ENDSECTION", "\n")
                yield KeyValue(line.strip(), "\n")

            # SUBJECT_SAMPLE_FACTORS line
            elif line.startswith("SUBJECT_SAMPLE_FACTORS"):
                line_items = line.split("\t")
                
                factor_pairs = line_items[3].split(" | ")
                factor_dict = dict_type()
                for pair in factor_pairs:
                    try:
                        factor_key, factor_value = pair.split(":")
                    except ValueError as e:
                        raise ValueError("Expected exactly 1 ':' in the factor key value pair, '" + pair + "'") from e
                    factor_key = factor_key.strip()
                    factor_value = factor_value.strip()
                    factor_dict[factor_key] = factor_value
                
                subject_sample_factors_dict = {
                    "Subject ID": line_items[1],
                    "Sample ID": line_items[2],
                    "Factors": factor_dict
                }
                if line_items[4]:
                    additional_data = dict_type()
                    for factor_item in line_items[4].split("; "):
                        try:
                            key, value = factor_item.split("=")
                        except ValueError as e:
                            raise ValueError("Expected exactly 1 '=' in the additional data key value pair, '" + factor_item + "'") from e
                        key = key.strip()
                        value = value.strip()
                        additional_data[key] = value
                    subject_sample_factors_dict["Additional sample data"] = additional_data
                yield KeyValue(line_items[0].strip(), subject_sample_factors_dict)

            # data start header
            elif line.endswith("_START"):
                yield KeyValue(line, "\n")

                # tokenize lines in data section till line ending with "_END" is reached
                while not line.endswith("_END"):
                    line = stream.popleft()
                    if line.endswith("_END"):
                        yield KeyValue(line.strip(), "\n")
                    else:
                        data = line.split("\t")
                        data = [token.strip('" ') for token in data]
                        yield KeyValue(data[0], tuple(data))

            # item line in item section (e.g. PROJECT, SUBJECT, etc..)
            elif line:
                if "_RESULTS_FILE" in line:
                    line_items = line.split("\t")
                    yield KeyValue(line_items[0].strip()[3:], _results_file_line_to_dict('\t'.join(line_items[1:])))
                else:
                    try:
                        key, value = line.split("\t", 1)
                    except ValueError as e:
                        raise ValueError("Expected a tab in the line.") from e
                    if ":" in key:
                        if ":UNITS" in key:
                            yield KeyValue("Units", value)
                        else:
                            yield KeyValue(key.strip()[3:], value)
                    else:
                        yield KeyValue(key.strip(), value.strip())

        except IndexError as e:
            raise IndexError("LINE WITH ERROR:\n\t" + repr(line)) from e
        except ValueError as e:
            raise ValueError("LINE WITH ERROR:\n\t" + repr(line)) from e

    # end of file
    yield KeyValue("#ENDSECTION", "\n")
    yield KeyValue("!#ENDFILE", "\n")  # This is to ensure that tokenizer terminates when #END is missing.
