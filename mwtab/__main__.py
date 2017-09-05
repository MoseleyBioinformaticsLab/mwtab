#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import docopt

from . import cli
from . import __version__


if __name__ == "__main__":

    args = docopt.docopt(cli.__doc__, version=__version__)
    cli.cli(args)
