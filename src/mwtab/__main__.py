#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import docopt

from . import cli
from . import __version__


def main():
    doc = [line for line in cli.__doc__.split('\n')]
    doc = doc[:3] + [line.lstrip() for line in doc[5:]]
    doc = '\n'.join(doc)
    args = docopt.docopt(cli.__doc__, version=__version__)
    cli.cli(args)


if __name__ == "__main__":
    main()
