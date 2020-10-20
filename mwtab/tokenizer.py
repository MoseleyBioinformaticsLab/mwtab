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


KeyValue = namedtuple("KeyValue", ["key", "value"])
SubjectSampleFactors = namedtuple("SubjectSampleFactors", ["key", "subject_type", "local_sample_id", "factors", "additional_sample_data"])
KeyValueExtra = namedtuple("KeyValueExtra", ["key", "value", "extra"])


def tokenizer(text, verbose=False):
    """A lexical analyzer for the `mwtab` formatted files.

    :param str text: `mwtab` formatted text.
    :return: Tuples of data.
    :rtype: py:class:`~collections.namedtuple`
    """

    stream = deque(text.split("\n"))

    while len(stream) > 0:
        line = stream.popleft()

        if line.startswith("#METABOLOMICS WORKBENCH"):
            yield KeyValue("#METABOLOMICS WORKBENCH", "\n")
            # yield KeyValue("HEADER", line)

            for identifier in line.split(" "):
                if ":" in identifier:
                    key, value = identifier.split(":")
                    yield KeyValue(key, value)

        elif line.startswith("#ANALYSIS TYPE"):
            yield KeyValue("HEADER", line)

        elif line.startswith("#SUBJECT_SAMPLE_FACTORS:"):
            yield KeyValue("#ENDSECTION", "\n")
            yield KeyValue("#SUBJECT_SAMPLE_FACTORS", "\n")

        elif line.startswith("#"):
            yield KeyValue("#ENDSECTION", "\n")
            yield KeyValue(line.strip(), "\n")

        elif line.startswith("SUBJECT_SAMPLE_FACTORS"):
            key, subject_type, local_sample_id, factors, additional_sample_data = line.split("\t")
            factors = {factor_item.split(":")[0].strip(): factor_item.split(":")[1].strip() for factor_item in factors.split("|")}
            additional_sample_dict = dict()
            # if additional_sample_data:
            #     additional_sample_data = {add_item.split("=")[0].strip(): add_item.split("=")[1].strip() for add_item in additional_sample_data.split(";")}
            for item in additional_sample_data.split(";"):
                if "=" in item:
                    key, value = item.split(":")
                    additional_sample_dict[key] = value
            yield SubjectSampleFactors(key.strip(), subject_type, local_sample_id, factors, additional_sample_dict)

        elif line.endswith("_START"):
            yield KeyValue(line, "\n")

            while not line.endswith("_END"):
                line = stream.popleft()
                if line.endswith("_END"):
                    yield KeyValue(line.strip(), "\n")
                else:
                    data = line.split("\t")
                    yield KeyValue(data[0], tuple(data))

        else:
            if line:
                if line.startswith("MS:MS_RESULTS_FILE") or line.startswith("NM:NMR_RESULTS_FILE"):
                    try:
                        key, value = line.split("\t")
                        yield KeyValue(key.strip()[3:], value)
                    except ValueError:
                        split_keys_values = line.split("\t")
                        key = split_keys_values[0]
                        value = split_keys_values[1]
                        extra = [tuple(pair.split(":")) for pair in split_keys_values[2:]]
                        yield KeyValueExtra(key.strip()[3:], value, extra)
                else:
                    try:
                        key, value = line.split("\t")
                        if ":" in key:
                            if key.startswith("MS_METABOLITE_DATA:UNITS"):
                                yield KeyValue(key.strip(), value)
                            else:
                                yield KeyValue(key.strip()[3:], value)
                        else:
                            yield KeyValue(key.strip(), value)
                    except ValueError:
                        if verbose:
                            print("LINE WITH ERROR:\n\t", repr(line))
                        raise ValueError("LINE WITH ERROR:\n\t", repr(line))

    yield KeyValue("#ENDSECTION", "\n")
    yield KeyValue("!#ENDFILE", "\n")  # This is to ensure that tokenizer terminates when #END is missing.
