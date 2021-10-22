#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Routines for working with ``mwTab`` format files used by the 
Metabolomics Workbench.

This package includes the following modules:

``mwtab``
    This module provides the :class:`~mwtab.mwtab.MWTabFile` class which is a python
    dictionary representation of a Metabolomics Workbench `mwtab` file. Data can be accessed
    directly from the :class:`~mwtab.mwtab.MWTabFile` instance using bracket accessors.

``cli``
    This module provides command-line interface for the ``mwtab`` package.

``tokenizer``
    This module provides the :func:`~mwtab.tokenizer.tokenizer` generator that generates
    tuples of key-value pairs from `mwtab` files.

``fileio``
    This module provides the :func:`~mwtab.fileio.read_files` generator
    to open files from different sources (single file/multiple files on a local 
    machine, directory/archive of files, URL address of a file).

``converter``
    This module provides the :class:`~mwtab.converter.Converter` class that is
    responsible for the conversion of ``mwTab`` formated files into their JSON
    representation and vice versa.

``mwschema``
    This module provides JSON schema definitions for the ``mwTab`` formatted files,
    i.e. specifies required and optional keys as well as data types.

``validator``
    This module provides routines to validate ``mwTab`` formatted files based
    on schema definitions as well as checks for file self-consistency.

``mwrest``
    This module provides the :class:`~mwtab.mwrest.GenericMWURL` class which is a
    python dictionary representation of a Metabolomics Workbench REST URL. The class
    is used to validate query parameters and to generate a URL path which can be
    used to request data from Metabolomics Workbench through their REST API.
"""

from logging import getLogger, NullHandler
from .fileio import read_files, read_mwrest
from .validator import validate_file
from .mwrest import GenericMWURL


__version__ = "1.2.2"


# Setting default logging handler
getLogger(__name__).addHandler(NullHandler())
