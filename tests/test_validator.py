import pytest
import mwtab


@pytest.mark.parametrize("files_source", [
    "tests/example_data/validation_files/AN000001.txt",
    "tests/example_data/validation_files/AN000041.txt"
])
def test_validate(files_source):
    mwfile = next(mwtab.read_files(files_source))
    validation_errors = mwtab.validate_file(mwfile)
    assert not validation_errors


def test_validate_ms_samples():
    mwfile = next(mwtab.read_files("example_data/validation_files/AN000001_error_1.txt"))
    validation_errors = mwtab.validate_file(mwfile, validate_factors=False, validate_features=False,
                                            validate_schema=False)
    assert len(validation_errors) == 1
    test_error = KeyError("Missing key `Samples` in `MS_METABOLITE_DATA` block.")
    assert type(validation_errors[0]) == type(test_error) and validation_errors[0].args == test_error.args

    mwfile = next(mwtab.read_files("example_data/validation_files/AN000001_error_2.txt"))
    validation_errors = mwtab.validate_file(mwfile, validate_factors=False, validate_features=False,
                                            validate_schema=False)
    assert len(validation_errors) == 3
    test_error = ValueError("Sample with no Sample ID (\"\") in `SUBJECT_SAMPLE_FACTOR` block.")
    assert type(validation_errors[0]) == type(test_error) and validation_errors[0].args == test_error.args
    test_error = ValueError("Sample with no Sample ID (\"\") in `MS_METABOLITE_DATA` block (usually caused by "
                            "extraneous TAB at the end of line).")
    assert type(validation_errors[1]) == type(test_error) and validation_errors[1].args == test_error.args
    test_error = ValueError("`MS_METABOLITE_DATA` block contains additional samples not found in "
                            "`SUBJECT_SAMPLE_FACTORS` block.\n\tAdditional samples: {'LabF_115873'}")
    assert type(validation_errors[2]) == type(test_error) and validation_errors[2].args == test_error.args


def test_validate_nmr_samples():
    mwfile = next(mwtab.read_files("example_data/validation_files/AN000041_error_1.txt"))
    validation_errors = mwtab.validate_file(mwfile, validate_factors=False, validate_features=False,
                                            validate_schema=False)
    assert len(validation_errors) == 1
    test_error = KeyError("Missing key `Bin range(ppm)` in `NMR_BINNED_DATA` block.")
    assert type(validation_errors[0]) == type(test_error) and validation_errors[0].args == test_error.args

    mwfile = next(mwtab.read_files("example_data/validation_files/AN000041_error_2.txt"))
    validation_errors = mwtab.validate_file(mwfile, validate_factors=False, validate_features=False,
                                            validate_schema=False)
    assert len(validation_errors) == 3
    test_error = ValueError("Sample with no Sample ID (\"\") in `SUBJECT_SAMPLE_FACTOR` block.")
    assert type(validation_errors[0]) == type(test_error) and validation_errors[0].args == test_error.args
    test_error = ValueError("Sample with no sample ID (\"\") in `NMR_BINNED_DATA` block (usually caused by extraneous "
                            "TAB at the end of line).")
    assert type(validation_errors[1]) == type(test_error) and validation_errors[1].args == test_error.args
    test_error = ValueError("`NMR_BINNED_DATA` block contains additional samples not found in `SUBJECT_SAMPLE_FACTORS` "
                            "block.\n\tAdditional samples: {'C0559'}")
    assert type(validation_errors[2]) == type(test_error) and validation_errors[2].args == test_error.args


def test_validate_factors():
    mwfile = next(mwtab.read_files("example_data/validation_files/AN000001_error_3.txt"))
    validation_errors = mwtab.validate_file(mwfile, validate_samples=False, validate_features=False,
                                            validate_schema=False)
    assert len(validation_errors) == 1
    assert repr(validation_errors[0]) == "KeyError('Missing key `Factors` in `MS_METABOLITE_DATA` block.')"

    # mwfile = next(mwtab.read_files("example_data/validation_files/AN000001_error_4.txt"))
    # validation_errors = mwtab.validate_file(mwfile, validate_samples=False, validate_features=False, validate_schema=False)
    # print(len(validation_errors))
    # print(validation_errors)
    # assert len(validation_errors) == 3
    # assert repr(validation_errors[0]) == "ValueError('Sample with no Factor(s) (\"\") in `SUBJECT_SAMPLE_FACTOR` block.',)"
    # assert repr(validation_errors[2]) == "ValueError('Sample with no factors (\"\") in `MS_METABOLITE_DATA` block " \
    #                                      "(usually caused by extraneous TAB at the end of line).',)"
    # assert repr(validation_errors[3]).replace("\"", "'").replace("\\n", "").replace("\\t", "") == "ValueError('`NMR_BINNED_DATA` block contains additional samples not found in `SUBJECT_SAMPLE_FACTORS` block.Additional samples: {'C0559'}',)"


def test_validate_metabolites():
    mwfile = next(mwtab.read_files("example_data/validation_files/AN000001_error_5.txt"))
    validation_errors = mwtab.validate_file(mwfile, validate_samples=False, validate_factors=False,
                                            validate_schema=False)
    assert len(validation_errors) == 1

    pass


processing_errors = [
    "AN000404", "AN000405", "AN000410", "AN000415", "AN000436", "AN000439", "AN001311", "AN001312", "AN001313",
    "AN000696", "AN001467", "AN001576", "AN001684", "AN001685", "AN001721", "AN001979", "AN002035", "AN001492",
    "AN001493", "AN001499", "AN001761", "AN001762", "AN001776", "AN001777", "AN001982", "AN001689", "AN001690",
    "AN001992"
]

REGEXS = [
    (r"(?i)(m/z)", "m/z"),                                          # m/z
    (r"(?i)(quan)[\S]{,}(\s|_)(m)[\S]{,}(z)", "quantitated_m/z"),   # quantitated_m/z
    (r"(?i)(r)[\s|\S]{,}(time)[\S]{,}", "retention_time"),          # retention_time
    (r"(?i)(ret)[\s|\S]{,}(index)", "retention_index"),             # retention_index
    (r"(?i)(pubchem)[\S]{,}", "pubchem_id"),                        # pubchem_id
    (r"(?i)(inchi)[\S]{,}", "inchi_key"),                           # inchi_key
    (r"(?i)(moverz)(\s|_)(quant)", "moverz_quant"),                 # moverz_quant
    (r"(?i)(kegg)(\s|_)(i)", "kegg_id"),                            # kegg_id
    (r"(?i)(kegg)$", "kegg_id"),
    (r"(?i)(ri)$", "ri"),                                           # ri
    (r"(?i)(ri)(\s|_)(type)", "ri_type"),                           # ri_type
    (r"(?i)(other)(\s|_)(id)", "other_id"),                         # other_id (other_id_type)
    (r"(?i)[\s|\S]{,}(HMDB)", "hmdb_id"),                           # hmdb_id
    (r"(?i)(Human Metabolome D)[\S]{,}", "hm0db_id"),
]

duplicate_fields = {f: dict() for r, f in REGEXS}

if __name__ == '__main__':

    test_validate("example_data/validation_files/AN000001.txt")
    test_validate("example_data/validation_files/AN000041.txt")
    test_validate_ms_samples()
    test_validate_nmr_samples()
    test_validate_factors()
    exit()

    # from os import walk
    # from re import match
    # from mwtab import read_files, se   ction_schema_mapping
    #
    # error_files = dict()
    # unique_fields = dict()
    # (_, _, filenames) = next(walk("/mlab/data/cdpo224/mwtab/data"))
    # filenames = sorted(filenames)
    # for filename in filenames:
    #     if not any(error in filename for error in processing_errors):
    #         mwfile = next(read_files("/mlab/data/cdpo224/mwtab/data/{}".format(filename)))
    #         if mwfile.get("METABOLITES"):
    #             if mwfile["METABOLITES"]["METABOLITES_START"].get("Fields"):
    #                 from_metabolites_fields = set(mwfile["METABOLITES"]["METABOLITES_START"]["Fields"])
    #                 for field in from_metabolites_fields:
    #                     if field in unique_fields.keys():
    #                         unique_fields[field] += 1
    #                     else:
    #                         unique_fields[field] = 1
    # del unique_fields["metabolite_name"]
    #
    # print(len(unique_fields))
    # sent = 0
    # items = list(unique_fields.items())
    # for k, v in items:
    #     for r, f in REGEXS:
    #         if match(r, k):
    #             duplicate_fields[f].update({k: v})
    #             del unique_fields[k]
    #             break
    #
    # features = sorted(unique_fields.items(), key=lambda x: x[0].lower(), reverse=True)
    # max_len = max(len(f[0]) for f in features)
    # print("Unique Features: ")
    # print("{0:<{1}} {2}".format("FEATURE", max_len, "INSTANCES"))
    # for f, n in features:
    #     print("{0:<{1}} {2}".format(f, max_len, n))
    #
    # print("\nDuplicate Features: ")
    # for k in duplicate_fields.keys():
    #     print(k)
    #     for f in duplicate_fields[k].keys():
    #         print("\t{}".format(f))
