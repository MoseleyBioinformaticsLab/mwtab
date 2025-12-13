#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
mwtab.fileio
~~~~~~~~~~~~

This module provides routines for reading ``mwTab`` formatted files
from difference kinds of sources:

   * Single ``mwTab`` formatted file on a local machine.
   * Directory containing multiple ``mwTab`` formatted files.
   * Compressed zip/tar archive of ``mwTab`` formatted files.
   * URL address of ``mwTab`` formatted file.
   * ``ANALYSIS_ID`` of ``mwTab`` formatted file. 
"""

import os
import io
import zipfile
import tarfile
import bz2
import gzip
from re import match
import pathlib
from typing import Any
from functools import partial

from . import mwtab
from . import mwrest

from urllib.request import urlopen
from urllib.parse import urlparse


VERBOSE = False
MWREST_URL = mwrest.BASE_URL


def _create_save_path(path):
    """Create directories in the path that don't already exist.
    
    :param str path: Path to save something to.
    """
    path = pathlib.Path(path)
    # Assume if there is a suffix (file extension) that it's a file and only make the parent path.
    if path.suffix:
        path.parent.mkdir(parents=True, exist_ok=True)
    else:
        path.mkdir(parents=True, exist_ok=True)


def _return_correct_yield(*payload, exception, return_exceptions):
    """Simple helper to yield either payload or payload and exception for generators.
    
    :param Any payload: Any object(s) to be yielded back.
    :param Exception or None exception: Any excpetion that occurred during the execution of the generator
    :param bool return_exceptions: Whether to yield a tuple with payload and exception or just the payload.
    :return: If return_exceptions is True return a tuple of payload and exception, else return only payload.
    """
    if exception is not None and not return_exceptions:
        raise exception
    return (*payload, exception) if return_exceptions else payload[0] if len(payload) == 1 else payload


def _generate_filenames(sources, return_exceptions=False):
    """Generate filenames.

    :param tuple sources: Sequence of strings representing path to file(s).
    :param bool return_exceptions: Whether to yield a tuple with path and exception or just the path.
    :return: Path to file(s), or a tuple of path and exceptions.
    :rtype: :py:class:`str`
    """
    for source in sources:
        try:
            if os.path.isdir(source):
                for path, _, filelist in os.walk(source):
                    for fname in sorted(filelist):
                        if os.path.splitext(fname)[1].lower() in {".csv", ".txt", ".json"}:
                            yield _return_correct_yield(os.path.join(path, fname), 
                                                        exception=None, 
                                                        return_exceptions=return_exceptions)
    
            elif os.path.isfile(source):
                yield _return_correct_yield(source, 
                                            exception=None, 
                                            return_exceptions=return_exceptions)
    
            elif source.isdigit():
                url, e = next(mwrest.generate_mwtab_urls([source], base_url=MWREST_URL, return_exceptions=True))
                yield _return_correct_yield(url, 
                                            exception=e, 
                                            return_exceptions=return_exceptions)
    
            # TODO: Add ST parsing
            elif match(r"(AN[0-9]{6}$)", source):
                url, e = next(mwrest.generate_mwtab_urls([source], base_url=MWREST_URL, return_exceptions=True))
                yield _return_correct_yield(url, 
                                            exception=e, 
                                            return_exceptions=return_exceptions)
    
            elif GenericFilePath.is_url(source):
                yield _return_correct_yield(source, 
                                            exception=None, 
                                            return_exceptions=return_exceptions)
    
            else:
                yield _return_correct_yield(source, 
                                            exception=TypeError("Unknown file source."), 
                                            return_exceptions=return_exceptions)
        except Exception as e:
            yield _return_correct_yield(source, 
                                        exception=e, 
                                        return_exceptions=return_exceptions)


def _generate_handles(filenames, return_exceptions=False):
    """Open a sequence of filenames one at time producing file objects.
    The file is closed immediately when proceeding to the next iteration.

    :param generator filenames: Generator object that yields the path to each file, one at a time.
    :param bool return_exceptions: Whether to yield a tuple with filehandle and exception or just the filehandle.
    :return: Filehandle to be processed into an instance.
    """
    for fname, exc in filenames:
        if exc is not None:
            yield _return_correct_yield(None, fname,
                                        exception=exc, 
                                        return_exceptions=return_exceptions)
            continue
        try:
            path = GenericFilePath(fname)
            for filehandle, source in path.open():
                yield _return_correct_yield(filehandle, source, 
                                            exception=None, 
                                            return_exceptions=return_exceptions)
                filehandle.close()
        except Exception as e:
            yield _return_correct_yield(None, fname, 
                                        exception=e, 
                                        return_exceptions=return_exceptions)


def read_with_class(sources: str|list[str], read_class: type, class_kwds: dict, return_exceptions: bool = False) -> tuple[Any, Exception]|Any:
    """Read from sources using the given read_class.
    
    This is really created to use functools partial to create a read mwthod for a particular class.
    
    Args:
        sources: A string or list of strings to read from.
        read_class: A class with a read() method to instantiate to read from source.
        class_kwds: A dictionary of keyword arguments to pass to the class constructor.
        return_exceptions: Whether to yield a tuple with file instance and exception or just the file instance.
    
    Returns:
        Returns the instantiated class and any exceptions, or None and any exceptions, or the source and any exceptions.
    """
    sources = [sources] if not isinstance(sources, list) else sources
    try:
        filenames = _generate_filenames(sources, True)
        filehandles = _generate_handles(filenames, True)
    except Exception as e:
        yield _return_correct_yield(None, 
                                    exception=e, 
                                    return_exceptions=return_exceptions)
    for fh, source, exc in filehandles:
        if exc is not None:
            yield _return_correct_yield(source, 
                                        exception=exc, 
                                        return_exceptions=return_exceptions)
            continue
        try:
            f = read_class(source, **class_kwds)
            f.read(fh)
            fh.close()

            if VERBOSE:
                print("Processed file: {}".format(os.path.abspath(source)))
            
            yield _return_correct_yield(f, 
                                        exception=None, 
                                        return_exceptions=return_exceptions)

        except Exception as e:
            fh.close()
            if VERBOSE:
                print("Error processing file: ", os.path.abspath(source), "\nReason:", e)
            yield _return_correct_yield(source, 
                                        exception=e, 
                                        return_exceptions=return_exceptions)


read_files = partial(read_with_class, read_class = mwtab.MWTabFile, class_kwds = {"duplicate_keys": True})
read_mwrest = partial(read_with_class, read_class = mwrest.MWRESTFile, class_kwds = {})

class ReadLines():
    def __init__(self, source, *args, **kwargs):
        self.source = source
        self.lines = []
    
    def read(self, filehandle):
        """
        """
        string = filehandle.read()
        filehandle.close()
        if isinstance(string, str):
            lines = string.replace("\r", "\n").split("\n")
        elif isinstance(string, bytes):
            lines = string.decode("utf-8").replace("\r", "\n").split("\n")
        else:
            raise TypeError("Expecting <class 'str'> or <class 'bytes'>, but {} was passed".format(type(string)))
        
        self.lines = [line for line in lines if line]
        

read_lines = partial(read_with_class, read_class = ReadLines, class_kwds = {})


class GenericFilePath(object):
    """`GenericFilePath` class knows how to open local files or files over URL."""

    def __init__(self, path):
        """Initialize path.

        :param str path: String representing a path to local file(s) or valid URL address of file(s).
        """
        self.path = path

    def open(self):
        """Generator that opens and yields filehandles using appropriate facilities:
        test if path represents a local file or file over URL, if file is compressed
        or not.

        :return: Filehandle to be processed into an instance.
        """
        is_url = self.is_url(self.path)
        compression_type = self.is_compressed(self.path)

        if not compression_type:
            if is_url:
                filehandle = urlopen(self.path)
            else:
                filehandle = open(self.path, "r", encoding="utf-8")
            source = self.path
            yield filehandle, source
            filehandle.close()

        elif compression_type:
            if is_url:
                response = urlopen(self.path)
                path = response.read()
                response.close()
            else:
                path = self.path

            if compression_type == "zip":
                ziparchive = zipfile.ZipFile(io.BytesIO(path), "r") if is_url else zipfile.ZipFile(path)
                for name in ziparchive.infolist():
                    if not name.filename.endswith("/"):
                        filehandle = ziparchive.open(name)
                        source = self.path + "/" + name.filename
                        yield filehandle, source
                        filehandle.close()

            elif compression_type in ("tar", "tar.bz2", "tar.gz"):
                tararchive = tarfile.open(fileobj=io.BytesIO(path)) if is_url else tarfile.open(path)
                for name in tararchive:
                    if name.isfile():
                        filehandle = tararchive.extractfile(name)
                        source = self.path + "/" + name.name
                        yield filehandle, source
                        filehandle.close()

            elif compression_type == "bz2":
                filehandle = bz2.BZ2File(io.BytesIO(path)) if is_url else bz2.BZ2File(path)
                source = self.path
                yield filehandle, source
                filehandle.close()

            elif compression_type == "gz":
                filehandle = gzip.open(io.BytesIO(path)) if is_url else gzip.open(path)
                source = self.path
                yield filehandle, source
                filehandle.close()

    @staticmethod
    def is_compressed(path):
        """Test if path represents compressed file(s).

        :param str path: Path to file(s).
        :return: String specifying compression type if compressed, "" otherwise.
        :rtype: :py:class:`str`
        """
        if path.endswith(".zip"):
            return "zip"
        elif path.endswith(".tar.gz"):
            return "tar.gz"
        elif path.endswith(".tar.bz2"):
            return "tar.bz2"
        elif path.endswith(".gz"):
            return "gz"
        elif path.endswith(".bz2"):
            return "bz2"
        elif path.endswith(".tar"):
            return "tar"
        return ""

    @staticmethod
    def is_url(path):
        """Test if path represents a valid URL.

        :param str path: Path to file.
        :return: True if path is valid url string, False otherwise.
        :rtype: :py:obj:`True` or :py:obj:`False`
        """
        try:
            parse_result = urlparse(path)
            return all((parse_result.scheme, parse_result.netloc, parse_result.path))
        except ValueError:
            return False
