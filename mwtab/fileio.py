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
import sys
import zipfile
import tarfile
import bz2
import gzip

from . import mwtab
from . import validator
from . import mwschema


if sys.version_info.major == 3:
    from urllib.request import urlopen
    from urllib.parse import urlparse
else:
    from urllib2 import urlopen
    from urlparse import urlparse


VERBOSE = False
MWREST = "http://www.metabolomicsworkbench.org/rest/study/analysis_id/{}/mwtab/txt"


def _generate_filenames(sources):
    """Generate filenames.

    :param tuple sources: Sequence of strings representing path to file(s).
    :return: Path to file(s).
    :rtype: :py:class:`str`
    """
    for source in sources:
        if os.path.isdir(source):
            for path, dirlist, filelist in os.walk(source):
                for fname in filelist:
                    if GenericFilePath.is_compressed(fname):
                        if VERBOSE:
                            print("Skipping compressed file: {}".format(os.path.abspath(fname)))
                        continue
                    else:
                        yield os.path.join(path, fname)

        elif os.path.isfile(source):
            yield source

        elif source.isdigit():
            analysis_id = "AN{}".format(source.zfill(6))
            url = MWREST.format(analysis_id)
            yield url

        elif GenericFilePath.is_url(source):
            yield source

        else:
            raise TypeError("Unknown file source.")


def _generate_handles(filenames):
    """Open a sequence of filenames one at time producing file objects.
    The file is closed immediately when proceeding to the next iteration.

    :param generator filenames: Generator object that yields the path to each file, one at a time.
    :return: Filehandle to be processed into an instance.
    """
    for fname in filenames:
        path = GenericFilePath(fname)
        for filehandle, source in path.open():
            yield filehandle, source
            filehandle.close()


def read_files(*sources, **kwds):
    """Construct a generator that yields file instances.

    :param sources: One or more strings representing path to file(s).
    """
    filenames = _generate_filenames(sources)
    filehandles = _generate_handles(filenames)
    for fh, source in filehandles:
        try:
            f = mwtab.MWTabFile(source)
            f.read(fh)

            if kwds.get('validate'):
                validator.validate_file(mwtabfile=f,
                                        section_schema_mapping=mwschema.section_schema_mapping,
                                        validate_samples=True,
                                        validate_factors=True)
            yield f

            if VERBOSE:
                print("Processed file: {}".format(os.path.abspath(source)))

        except Exception as e:
            if VERBOSE:
                print("Error processing file: ", os.path.abspath(source), "\nReason:", e)
            pass


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
                filehandle = open(self.path, "r")
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
