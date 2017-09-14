#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
mwtab.converter
~~~~~~~~~~~~~~~

This module provides functionality for converting between the
Metabolomics Workbench ``mwTab`` formatted file and its equivalent
JSONized representation.

The following conversions are possible:

Local files:
   * One-to-one file conversions:
      * textfile - to - textfile
      * textfile - to - textfile.gz
      * textfile - to - textfile.bz2
      * textfile.gz - to - textfile
      * textfile.gz - to - textfile.gz
      * textfile.gz - to - textfile.bz2
      * textfile.bz2 - to - textfile
      * textfile.bz2 - to - textfile.gz
      * textfile.bz2 - to - textfile.bz2
      * textfile / textfile.gz / textfile.bz2 - to - textfile.zip / textfile.tar / textfile.tar.gz / textfile.tar.bz2 (TypeError: One-to-many conversion)
   * Many-to-many files conversions:
      * Directories:
         * directory - to - directory
         * directory - to - directory.zip
         * directory - to - directory.tar
         * directory - to - directory.tar.bz2
         * directory - to - directory.tar.gz
         * directory - to - directory.gz / directory.bz2 (TypeError: Many-to-one conversion)
      * Zipfiles:
         * zipfile.zip - to - directory
         * zipfile.zip - to - zipfile.zip
         * zipfile.zip - to - tarfile.tar
         * zipfile.zip - to - tarfile.tar.gz
         * zipfile.zip - to - tarfile.tar.bz2
         * zipfile.zip - to - directory.gz / directory.bz2 (TypeError: Many-to-one conversion)
      * Tarfiles:
         * tarfile.tar - to - directory
         * tarfile.tar - to - zipfile.zip
         * tarfile.tar - to - tarfile.tar
         * tarfile.tar - to - tarfile.tar.gz
         * tarfile.tar - to - tarfile.tar.bz2
         * tarfile.tar - to - directory.gz / directory.bz2 (TypeError: Many-to-one conversion)
         * tarfile.tar.gz - to - directory
         * tarfile.tar.gz - to - zipfile.zip
         * tarfile.tar.gz - to - tarfile.tar
         * tarfile.tar.gz - to - tarfile.tar.gz
         * tarfile.tar.gz - to - tarfile.tar.bz2
         * tarfile.tar.gz - to - directory.gz / directory.bz2 (TypeError: Many-to-one conversion)
         * tarfile.tar.bz2 - to - directory
         * tarfile.tar.bz2 - to - zipfile.zip
         * tarfile.tar.bz2 - to - tarfile.tar
         * tarfile.tar.bz2 - to - tarfile.tar.gz
         * tarfile.tar.bz2 - to - tarfile.tar.bz2
         * tarfile.tar.bz2 - to - directory.gz / directory.bz2 (TypeError: Many-to-one conversion)
URL files:
   * One-to-one file conversions:
      * analysis_id - to - textfile
      * analysis_id - to - textfile.gz
      * analysis_id - to - textfile.bz2
      * analysis_id - to - textfile.zip / textfile.tar / textfile.tar.gz / textfile.tar.bz2 (TypeError: One-to-many conversion)
      * textfileurl - to - textfile
      * textfileurl - to - textfile.gz
      * textfileurl - to - textfile.bz2
      * textfileurl.gz - to - textfile
      * textfileurl.gz - to - textfile.gz
      * textfileurl.gz - to - textfile.bz2
      * textfileurl.bz2 - to - textfile
      * textfileurl.bz2 - to - textfile.gz
      * textfileurl.bz2 - to - textfile.bz2
      * textfileurl / textfileurl.gz / textfileurl.bz2 - to - textfile.zip / textfile.tar / textfile.tar.gz / textfile.tar.bz2 (TypeError: One-to-many conversion)
   * Many-to-many files conversions:
      * Zipfiles:
         * zipfileurl.zip - to - directory
         * zipfileurl.zip - to - zipfile.zip
         * zipfileurl.zip - to - tarfile.tar
         * zipfileurl.zip - to - tarfile.tar.gz
         * zipfileurl.zip - to - tarfile.tar.bz2
         * zipfileurl.zip - to - directory.gz / directory.bz2 (TypeError: Many-to-one conversion)
      * Tarfiles:
         * tarfileurl.tar - to - directory
         * tarfileurl.tar - to - zipfile.zip
         * tarfileurl.tar - to - tarfile.tar
         * tarfileurl.tar - to - tarfile.tar.gz
         * tarfileurl.tar - to - tarfile.tar.bz2
         * tarfileurl.tar - to - directory.gz / directory.bz2 (TypeError: Many-to-one conversion)
         * tarfileurl.tar.gz - to - directory
         * tarfileurl.tar.gz - to - zipfile.zip
         * tarfileurl.tar.gz - to - tarfile.tar
         * tarfileurl.tar.gz - to - tarfile.tar.gz
         * tarfileurl.tar.gz - to - tarfile.tar.bz2
         * tarfileurl.tar.gz - to - directory.gz / directory.bz2 (TypeError: Many-to-one conversion)
         * tarfileurl.tar.bz2 - to - directory
         * tarfileurl.tar.bz2 - to - zipfile.zip
         * tarfileurl.tar.bz2 - to - tarfile.tar
         * tarfileurl.tar.bz2 - to - tarfile.tar.gz
         * tarfileurl.tar.bz2 - to - tarfile.tar.bz2
         * tarfileurl.tar.bz2 - to - directory.gz / directory.bz2 (TypeError: Many-to-one conversion)
"""

import os
import io
import zipfile
import tarfile
import bz2
import gzip

from . import fileio


class Translator(object):
    """Translator abstract class."""

    def __init__(self, from_path, to_path, from_format=None, to_format=None, validate=False):
        """Translator initializer.
        :param str from_path: Path to input file(s).
        :param str to_path: Path to output file(s).
        :param str from_format: Input format.
        :param str to_format: Output format.
        """
        self.from_path = from_path
        self.to_path = to_path
        self.from_format = from_format
        self.to_format = to_format
        self.from_path_compression = fileio.GenericFilePath.is_compressed(from_path)
        self.to_path_compression = fileio.GenericFilePath.is_compressed(to_path)
        self.validate = validate

    def __iter__(self):
        """Abstract iterator must be implemented in a subclass."""
        raise NotImplementedError()


class MWTabFileToMWTabFile(Translator):
    """Translator concrete class that can convert between ``mwTab`` and ``JSON`` formats."""

    file_extension = {"json": ".json",
                      "mwtab": ".txt"}

    def __init__(self, from_path, to_path, from_format=None, to_format=None, validate=False):
        """StarFileToStarFile translator initializer.
        :param str from_path: Path to input file(s).
        :param str to_path: Path to output file(s).
        :param str from_format: Input format: `mwtab` or `json`.
        :param str to_format: Output format: `mwtab` or `json`.
        """
        super(MWTabFileToMWTabFile, self).__init__(from_path, to_path, from_format, to_format, validate)

    def __iter__(self):
        """Iterator that yields instances of :class:`~mwtab.mwtab.MWTabFile` instances.
        :return: instance of :class:`~mwtab.mwtab.MWTabFile` object instance.
        :rtype: :class:`~mwtab.mwtab.MWTabFile`
        """
        for mwtabfile in fileio.read_files(self.from_path, validate=self.validate):
            yield mwtabfile


class Converter(object):
    """Converter class to convert ``mwTab`` files from ``mwTab`` to ``JSON`` or from ``JSON`` to ``mwTab`` format."""

    def __init__(self, from_path, to_path, from_format="mwtab", to_format="json", validate=False):
        """Converter initializer.
        :param file_generator:
        :type file_generator: :class:`mwtab.converter.Translator`
        """
        self.file_generator = MWTabFileToMWTabFile(from_path, to_path, from_format, to_format, validate)

    def convert(self):
        """Convert file(s) from ``mwTab`` format to ``JSON`` format or from ``JSON`` format to ``mwTab`` format.
        :return: None
        :rtype: :py:obj:`None`
        """
        if not os.path.exists(os.path.dirname(self.file_generator.to_path)):
            dirname = os.path.dirname(self.file_generator.to_path)
            if dirname:
                os.makedirs(dirname)

        if os.path.isdir(self.file_generator.from_path):
            self._many_to_many()
        elif os.path.isfile(self.file_generator.from_path) or fileio.GenericFilePath.is_url(self.file_generator.from_path):
            if self.file_generator.from_path_compression in ("zip", "tar", "tar.gz", "tar.bz2"):
                self._many_to_many()
            elif self.file_generator.from_path_compression in ("gz", "bz2"):
                self._one_to_one()
            elif not self.file_generator.from_path_compression:
                self._one_to_one()
        elif self.file_generator.from_path.isdigit():
            self._one_to_one()
        else:
            raise TypeError('Unknown input file format: "{}"'.format(self.file_generator.from_path))

    def _many_to_many(self):
        """Perform many-to-many files conversion.
        :return: None
        :rtype: :py:obj:`None`
        """
        if not self.file_generator.to_path_compression:
            self._to_dir(self.file_generator)
        elif self.file_generator.to_path_compression == "zip":
            self._to_zipfile(self.file_generator)
        elif self.file_generator.to_path_compression in ("tar", "tar.gz", "tar.bz2"):
            self._to_tarfile(self.file_generator)
        elif self.file_generator.to_path_compression in ("gz", "bz2"):
            raise TypeError('Many-to-one conversion, cannot convert "{}" into "{}"'.format(self.file_generator.from_path,
                                                                                           self.file_generator.to_path))
        else:
            raise TypeError('Unknown output file format: "{}"'.format(self.file_generator.to_path))

    def _one_to_one(self):
        """Perform one-to-one file conversion.
        :return: None
        :rtype: :py:obj:`None`
        """
        if not self.file_generator.to_path_compression:
            self._to_textfile(self.file_generator)
        elif self.file_generator.to_path_compression == "gz":
            self._to_gzipfile(self.file_generator)
        elif self.file_generator.to_path_compression == "bz2":
            self._to_bz2file(self.file_generator)
        elif self.file_generator.to_path_compression in ("tar", "tar.gz", "tar.bz2", "zip"):
            raise TypeError('One-to-many conversion, cannot convert "{}" into "{}"'.format(self.file_generator.from_path,
                                                                                           self.file_generator.to_path))
        else:
            raise TypeError('Unknown format: "{}"'.format(self.file_generator.to_path))

    def _to_dir(self, file_generator):
        """Convert files to directory.
        :return: None
        :rtype: :py:obj:`None`
        """
        for f in file_generator:
            outpath = self._output_path(f.source, file_generator.to_format)

            if not os.path.exists(os.path.dirname(outpath)):
                os.makedirs(os.path.dirname(outpath))

            with open(outpath, mode="w") as outfile:
                f.write(outfile, file_generator.to_format)

    def _to_zipfile(self, file_generator):
        """Convert files to zip archive.
        :return: None
        :rtype: :py:obj:`None`
        """
        with zipfile.ZipFile(file_generator.to_path, mode="w", compression=zipfile.ZIP_DEFLATED) as outfile:
            for f in file_generator:
                outpath = self._output_path(f.source, file_generator.to_format, archive=True)
                outfile.writestr(outpath, f.writestr(file_generator.to_format))

    def _to_tarfile(self, file_generator):
        """Convert files to tar archive.
        :return: None
        :rtype: :py:obj:`None`
        """
        if file_generator.to_path_compression == "tar":
            tar_mode = "w"
        elif file_generator.to_path_compression == "tar.gz":
            tar_mode = "w:gz"
        elif file_generator.to_path_compression == 'tar.bz2':
            tar_mode = "w:bz2"
        else:
            tar_mode = "w"

        with tarfile.open(file_generator.to_path, mode=tar_mode) as outfile:
            for f in file_generator:
                outpath = self._output_path(f.source, file_generator.to_format, archive=True)
                info = tarfile.TarInfo(outpath)
                data = f.writestr(file_generator.to_format).encode()
                info.size = len(data)
                outfile.addfile(tarinfo=info, fileobj=io.BytesIO(data))

    def _to_bz2file(self, file_generator):
        """Convert file to bz2-compressed file.
        :return: None
        :rtype: :py:obj:`None`
        """
        with bz2.BZ2File(file_generator.to_path, mode="wb") as outfile:
            for f in file_generator:
                outfile.write(f.writestr(file_generator.to_format).encode())

    def _to_gzipfile(self, file_generator):
        """Convert file to gzip-compressed file.
        :return: None
        :rtype: :py:obj:`None`
        """
        with gzip.GzipFile(file_generator.to_path, mode="wb") as outfile:
            for f in file_generator:
                outfile.write(f.writestr(file_generator.to_format).encode())

    def _to_textfile(self, file_generator):
        """Convert file to regular text file.
        :return: None
        :rtype: :py:obj:`None`
        """
        to_path = file_generator.to_path \
            if file_generator.to_path.endswith(file_generator.file_extension[file_generator.to_format]) \
            else file_generator.to_path + file_generator.file_extension[file_generator.to_format]

        with open(to_path, mode="w") as outfile:
            for f in file_generator:
                outfile.write(f.writestr(file_generator.to_format))

    def _output_path(self, input_path, to_format, archive=False):
        """Construct an output path string from an input path string.
        :param str input_path: Input path string.
        :return: Output path string.
        :rtype: :py:class:`str`
        """
        indirpath, fname = os.path.split(os.path.abspath(os.path.normpath(input_path)))

        commonprefix = os.path.commonprefix([os.path.abspath(self.file_generator.from_path),
                                             os.path.abspath(indirpath)])

        commonparts = commonprefix.split(os.sep)
        inparts = indirpath.split(os.sep)
        outparts = inparts[len(commonparts):]

        if archive:
            outdirpath = os.path.join(*outparts) if outparts else ""
        else:
            outdirpath = os.path.join(self.file_generator.to_path, *outparts)

        return os.path.join(outdirpath, fname + self.file_generator.file_extension[to_format])