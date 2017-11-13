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
"""

from .fileio import read_files

__version__ = "0.1.6"
