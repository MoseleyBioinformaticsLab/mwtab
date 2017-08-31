#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
mwtab.tokenizer
~~~~~~~~~~~~~~~

This module provides :func:`~mwtab.tokenizer.tokenizer` lexical analyzer for
`mwTab` format syntax. It is implemented as Python generator-based state
machine which generates (yields) tokens one at a time when :py:func:`next()`
is invoked on :func:`~mwtab.tokenizer.tokenizer` instance.

Each token is a tuple of "key-value"-like pairs, tuple of
``SUBJECT_SAMPLE_FACTORS`` or tuple of data deposited between
``*_START`` and ``*_END`` blocks.
"""

from collections import deque
from collections import namedtuple


KeyValue = namedtuple("KeyValue", ["key", "value"])
SubjectSampleFactors = namedtuple("SubjectSampleFactors", ["key", "subject_type", "local_sample_id", "factors", "additional_sample_data"])


def tokenizer(text):
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
            yield KeyValue("HEADER", line)

            for identifier in line.split(" "):
                if ":" in identifier:
                    key, value = identifier.split(":")
                    yield KeyValue(key, value)

        elif line.startswith("#SUBJECT_SAMPLE_FACTORS:"):
            yield KeyValue("#ENDSECTION", "\n")
            yield KeyValue("#SUBJECT_SAMPLE_FACTORS", "\n")

        elif line.startswith("#"):
            yield KeyValue("#ENDSECTION", "\n")
            yield KeyValue(line.strip(), "\n")

        elif line.startswith("SUBJECT_SAMPLE_FACTORS"):
            key, subject_type, local_sample_id, factors, additional_sample_data = line.split("\t")
            # factors = [dict([[i.strip() for i in f.split(":")]]) for f in factors.split("|")]
            yield SubjectSampleFactors(key.strip(), subject_type, local_sample_id, factors, additional_sample_data)

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
                try:
                    key, value, = line.split("\t")
                    if ":" in key:
                        yield KeyValue(key.strip()[3:], value)
                    else:
                        yield KeyValue(key.strip(), value)
                except Exception:
                    raise

    yield KeyValue("#ENDSECTION", "\n")
    yield KeyValue("!#ENDFILE", "\n")


if __name__ == "__main__":
    import sys

    script = sys.argv.pop(0)
    filepath = sys.argv.pop(0)

    with open(filepath, "r") as inf:
        t = tokenizer(inf.read())

    for token in t:
        print(token)
        pass
