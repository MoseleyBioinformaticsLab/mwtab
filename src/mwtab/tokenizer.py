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


KeyValue = namedtuple("KeyValue", ["key", "value"])
KeyValueExtra = namedtuple("KeyValueExtra", ["key", "value", "extra"])


def tokenizer(text):
    """A lexical analyzer for the `mwtab` formatted files.

    :param text: `mwTab` formatted text.
    :type text: py:class:`str`
    :return: Tuples of data.
    :rtype: py:class:`~collections.namedtuple`
    """
    stream = deque(text.split("\n"))

    while len(stream) > 0:
        line = stream.popleft()
        try:

            # header
            if line.startswith("#METABOLOMICS WORKBENCH"):
                yield KeyValue("#METABOLOMICS WORKBENCH", "\n")
                for identifier in line.split(" "):
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
                subject_sample_factors_dict = OrderedDict({
                    "Subject ID": line_items[1],
                    "Sample ID": line_items[2],
                    "Factors": {factor_item.split(":")[0].strip(): factor_item.split(":")[1].strip() for factor_item in
                                line_items[3].split("|")}
                })
                if line_items[4]:
                    subject_sample_factors_dict["Additional sample data"] = {
                        factor_item.split("=")[0].strip(): factor_item.split("=")[1].strip() for factor_item in line_items[4].split(";")
                    }
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
                        yield KeyValue(data[0], tuple(data))

            # item line in item section (e.g. PROJECT, SUBJECT, etc..)
            elif line:
                if "_RESULTS_FILE" in line:
                    line_items = line.split("\t")
                    # if len(line_items) > 2:
                    #     extra_items = list()
                    #     for extra_item in line_items[2:]:
                    #         k, v = extra_item.split(":")
                    #         extra_items.append(tuple([k.strip(), v.strip()]))
                    #     yield KeyValueExtra(line_items[0].strip()[3:], line_items[1], extra_items)
                    # else:
                    #     yield KeyValue(line_items[0].strip()[3:], line_items[1])
                    yield KeyValue(line_items[0].strip()[3:], " ".join(line_items[1:]))
                else:
                    key, value = line.split("\t")
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
