#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import json
from pprint import pprint
import docopt


from . import cli
from . import __version__
from . import mwtab
from . import mwschema
from . import validator

if __name__ == "__main__":
    # python3 -m mwtab tests/example_data/mwtab_studies/mwtab_files/ST000001_AN000001.txt

    script = sys.argv.pop(0)
    fpath = sys.argv.pop(0)

    with open(fpath, "r") as inf:
        mwtabfile = mwtab.MWTabFile(fpath)
        mwtabfile.read(inf)

    # mwtabfile.print_mwtabfile(file_format="mwtab", f=open("mf.txt", "w"))
    mwtabfile.print_file(file_format="mwtab")
    print("source:", mwtabfile.source)
    print("header:", mwtabfile.header)
    print("study id:", mwtabfile.study_id)
    print("analysis id:", mwtabfile.analysis_id)
    mwtabfile.print_file(file_format="json")

    # print(validator.validate_section(mwtabfile, "METABOLOMICS WORKBENCH", mwschema.metabolomics_workbench_schema))
    # print("\n")
    # print(validator.validate_section(mwtabfile, "PROJECT", mwschema.project_schema))
    # print("\n")
    # print(validator.validate_section(mwtabfile, "STUDY", mwschema.study_schema))
    # print("\n")
    # print(validator.validate_section(mwtabfile, "ANALYSIS", mwschema.analysis_schema))
    # print("\n")
    # print(validator.validate_section(mwtabfile, "SUBJECT_SAMPLE_FACTORS", mwschema.subject_sample_factors_schema))
    # print("\n")
    # print(validator.validate_section(mwtabfile, "COLLECTION", mwschema.collection_schema))
    # print("\n")
    # print(validator.validate_section(mwtabfile, "TREATMENT", mwschema.treatment_schema))
    # print("\n")
    # print(validator.validate_section(mwtabfile, "SAMPLEPREP", mwschema.sampleprep_schema))
    # print("\n")
    # print(validator.validate_section(mwtabfile, "MS", mwschema.ms_schema))
    # print("\n")
    # print(validator.validate_section(mwtabfile, "NMR", mwschema.nmr_schema))
    # print("\n")
    # print(validator.validate_section(mwtabfile, "METABOLITES", mwschema.metabolites_schema))
    # print("\n")
    # print(validator.validate_section(mwtabfile, "MS_METABOLITE_DATA", mwschema.ms_metabolite_data_schema))
    # print("\n")
    # print(validator.validate_section(mwtabfile, "NMR_BINNED_DATA", mwschema.nmr_binned_data_schema))

    try:
        validator.validate_file(mwtabfile=mwtabfile, section_schema_mapping=mwschema.section_schema_mapping)
    except Exception:
        raise

    #     # with open("{}.json".format(fpath), "w") as outf:
    #     #     json.dump(mwtabfile, outf, indent=4)
    #
    #     pprint(mwtabfile, width=200)

    # args = docopt.docopt(cli.__doc__, version=__version__)
    # cli.cli(args)
