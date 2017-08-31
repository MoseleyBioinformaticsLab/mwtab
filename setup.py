#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import re
from setuptools import setup, find_packages


if sys.argv[-1] == 'publish':
    os.system('python3 setup.py sdist upload')
    sys.exit()


def readme():
    with open('README.rst') as readme_file:
        return readme_file.read()


def find_version():
    with open('mwtab/__init__.py', 'r') as fd:
        version = re.search(r'^__version__\s*=\s*[\'"]([^\'"]*)[\'"]',
                            fd.read(), re.MULTILINE).group(1)
    if not version:
        raise RuntimeError('Cannot find version information')
    return version


REQUIRES = [
    "docopt >= 0.6.2",
    "schema >= 0.6.6"
]


setup(
    name='mwtab',
    version=find_version(),
    author='Andrey Smelter',
    author_email='andrey.smelter@gmail.com',
    description='Parser for mwtab files from the Metabolomics Workbench',
    keywords='mwtab metabolomics workbench',
    license='BSD',
    url='https://github.com/MoseleyBioinformaticsLab/mwtab',
    packages=find_packages(),
    platforms='any',
    long_description=readme(),
    install_requires=REQUIRES,
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ]
)