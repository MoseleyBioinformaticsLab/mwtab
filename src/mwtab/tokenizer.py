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
from collections import deque, namedtuple, OrderedDict
import re

import json_duplicate_keys as jdks

from .duplicates_dict import DuplicatesDict


KeyValue = namedtuple("KeyValue", ["key", "value"])
KeyValueExtra = namedtuple("KeyValueExtra", ["key", "value", "extra"])


def tokenizer(text, dict_type = None):
    """A lexical analyzer for the `mwtab` formatted files.

    :param text: `mwTab` formatted text.
    :type text: py:class:`str`
    :param dict_type: the type of dictionary to use, default is OrderedDict.
    :return: Tuples of data.
    :rtype: py:class:`~collections.namedtuple`
    """
    if dict_type is None:
        dict_type = OrderedDict
        
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
                # if compatability_mode:
                #     factor_dict = jdks.JSON_DUPLICATE_KEYS(OrderedDict())
                # else:
                #     factor_dict = OrderedDict()
                for pair in factor_pairs:
                    factor_key, factor_value = pair.split(":")
                    factor_key = factor_key.strip()
                    factor_value = factor_value.strip()
                    # if compatability_mode:
                    #     factor_dict.set(factor_key, factor_value, ordered_dict=True)
                    # else:
                    #     factor_dict[factor_key] = factor_value
                    factor_dict[factor_key] = factor_value
                
                subject_sample_factors_dict = OrderedDict({
                    "Subject ID": line_items[1],
                    "Sample ID": line_items[2],
                    "Factors": factor_dict
                })
                if line_items[4]:
                    additional_data = dict_type()
                    # if compatability_mode:
                    #     additional_data = jdks.JSON_DUPLICATE_KEYS(OrderedDict())
                    # else:
                    #     additional_data = OrderedDict()
                    for factor_item in line_items[4].split("; "):
                        key, value = factor_item.split("=")
                        key = key.strip()
                        value = value.strip()
                        # if compatability_mode:
                        #     additional_data.set(key, value, ordered_dict=True)
                        # else:
                        #     additional_data[key] = value
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
                    yield KeyValue(line_items[0].strip()[3:], line_items[1:])
                else:
                    key, value = line.split("\t", 1)
                    if ":" in key:
                        if ":UNITS" in key:
                            yield KeyValue("Units", value)
                        else:
                            yield KeyValue(key.strip()[3:], value)
                    else:
                        yield KeyValue(key.strip(), value)

        except IndexError as e:
            raise IndexError("LINE WITH ERROR:\n\t", repr(line), e)
        except ValueError as e:
            raise ValueError("LINE WITH ERROR:\n\t", repr(line), e)

    # end of file
    yield KeyValue("#ENDSECTION", "\n")
    yield KeyValue("!#ENDFILE", "\n")  # This is to ensure that tokenizer terminates when #END is missing.
