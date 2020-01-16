# TODO: This file should not be included in the public release
# TODO: update "if any(val for val in from_metabolite_data_features if val != ""):" lines to use diff.
import mwtab
from collections import OrderedDict
from os import walk
from re import match

"""
Error List
0 - Missing Sample ID in `SUBJECT_SAMPLE_FACTORS` block ("" entered as value).
1 - Missing Sample Factor(s) in `SUBJECT_SAMPLE_FACTORS` block ("" entered as value).

2 - Extra tab in `MS_METABOLITE_DATA` block (Samples line).
3 - Extra Sample ID in `MS_METABOLITE_DATA` block (Samples line).
4 - KeyError in `MS_METABOLITE_DATA` block (missing `Samples` key).

5 - Extra tab in `NMR_BINNED_DATA` block (Samples line).
6 - Extra Sample ID in `NMR_BINNED_DATA` block (Samples line).
7 - KeyError in `NMR_BINNED_DATA` block (missing `Fields` key).

8 - Extra tab in `MS_METABOLITE_DATA` block (Factors line).
9 - Extra Sample Factor in `MS_METABOLITE_DATA` block (Factors line).
10 - KeyError in `MS_METABOLITE_DATA` block (missing `Factors` key).

11 - Unequal `Features` in `METABOLITES` and `MS_METABOLITE_DATA` (features listed under `metabolite_name`).

14 - Extra keys in `METABOLITE` block (under `METABOLITES_START`, `DATA` keys).
15 - Missing key `metabolite_name` in `METABOLITE` block (commonly replaced with `Compound`).
"""

processing_errors = [
    "AN000404", "AN000405", "AN000410", "AN000415", "AN000436", "AN000439", "AN001311", "AN001312", "AN001313",
    "AN000696", "AN001467", "AN001576", "AN001684", "AN001685", "AN001721", "AN001979", "AN002035", "AN001492",
    "AN001493", "AN001499", "AN001761", "AN001762", "AN001776", "AN001777", "AN001982", "AN001689", "AN001690",
    "AN001992"
]


def _validate_samples_factors(mwtabfile, validate_samples=True, validate_factors=True):
    """Validate ``Samples`` and ``Factors`` identifiers across the file.

    :param mwtabfile: Instance of :class:`~mwtab.mwtab.MWTabFile`.
    :type mwtabfile: :class:`~mwtab.mwtab.MWTabFile`
    :param validate_samples:
    :type validate_samples: :py:obj:`bool`
    :param validate_factors:
    :type validate_factors: :py:obj:`bool`
    :return: List of errors (["None"] if no errors).
    :rtype: :py:obj:`list`
    """
    errors = []

    from_subject_samples = {i["local_sample_id"] for i in mwtabfile["SUBJECT_SAMPLE_FACTORS"]["SUBJECT_SAMPLE_FACTORS"]}
    from_subject_factors = {i["factors"] for i in mwtabfile["SUBJECT_SAMPLE_FACTORS"]["SUBJECT_SAMPLE_FACTORS"]}

    if "" in from_subject_samples:
        errors.append(ValueError("Sample with no sample ID (\"\") in `SUBJECT_SAMPLE_FACTOR` block."))
    if "" in from_subject_factors:
        errors.append(ValueError("Sample with no Factor(s) (\"\") in `SUBJECT_SAMPLE_FACTOR` block."))

    if validate_samples:

        if "MS_METABOLITE_DATA" in mwtabfile:
            try:
                from_metabolite_data_samples = set(mwtabfile["MS_METABOLITE_DATA"]["MS_METABOLITE_DATA_START"]["Samples"])
                if not from_metabolite_data_samples.issubset(from_subject_samples):
                    if "" in from_metabolite_data_samples:
                        errors.append(ValueError(
                            "Sample with no sample ID (\"\") in `MS_METABOLITE_DATA` block (usually caused by "
                            "excetrainious TAB at the end of line)."))
                    if any(val for val in from_metabolite_data_samples if val != ""):
                        errors.append(ValueError(
                            "`MS_METABOLITE_DATA` block contains additional samples not found in "
                            "`SUBJECT_SAMPLE_FACTORS` block.\n\t Additional samples: {}".format(
                                from_metabolite_data_samples.difference(from_subject_samples))))
            except KeyError:
                errors.append(KeyError("Missing key `Samples` in `MS_METABOLITE_DATA` block."))

        if "NMR_BINNED_DATA" in mwtabfile:
            try:
                from_nmr_binned_data_samples = set(mwtabfile["NMR_BINNED_DATA"]["NMR_BINNED_DATA_START"]["Fields"][1:])
                if not from_nmr_binned_data_samples.issubset(from_subject_samples):
                    if "" in from_nmr_binned_data_samples:
                        errors.append(ValueError(
                            "Sample with no sample ID (\"\") in `NMR_BINNED_DATA` block (usually caused by "
                            "excetrainious TAB at the end of line)."))
                    if any(val for val in from_nmr_binned_data_samples if val != ""):
                        errors.append(ValueError(
                            "`NMR_BINNED_DATA` block contains additional samples not found in `SUBJECT_SAMPLE_FACTORS` "
                            "block.\n\t Additional samples: {}".format(
                                from_nmr_binned_data_samples.difference(from_subject_samples))))
            except KeyError:
                errors.append(KeyError("Missing key `Samples` in `MS_METABOLITE_DATA` block."))

    if validate_factors:

        if "MS_METABOLITE_DATA" in mwtabfile:
            try:
                from_metabolite_data_factors = set(mwtabfile["MS_METABOLITE_DATA"]["MS_METABOLITE_DATA_START"]["Factors"])
                if not from_metabolite_data_factors.issubset(from_subject_factors):
                    if "" in from_metabolite_data_factors:
                        errors.append(ValueError(
                            "Sample with no factors (\"\") in `MS_METABOLITE_DATA` block (usually caused by "
                            "excetrainious TAB at the end of line)."))
                    if any(val for val in from_metabolite_data_factors if val != ""):
                        errors.append(ValueError(
                            "`MS_METABOLITE_DATA` block contains additional factors not found in "
                            "`SUBJECT_SAMPLE_FACTORS` block.\n\t Additional factors: {}".format(
                                from_metabolite_data_factors.difference(from_subject_samples))))
            except KeyError:
                errors.append(KeyError("Missing key `Factors` in `MS_METABOLITE_DATA` block."))

    return errors


def _validate_metabolites(mwtabfile, validate_features=True):
    """Validate metabolite ``Features`` across the file.

    :param mwtabfile: Instance of :class:`~mwtab.mwtab.MWTabFile`.
    :type mwtabfile: :class:`~mwtab.mwtab.MWTabFile`
    :param validate_features:
    :type validate_features: :py:obj:`bool`
    :return: List of errors (["None"] if no errors).
    :rtype: :py:obj:`list`
    """
    errors = []

    from_metabolites_features = {i["metabolite_name"] for i in mwtabfile["METABOLITES"]["METABOLITES_START"]["DATA"]}

    if "" in from_metabolites_features:
        errors.append(ValueError(
            "Feature with no name (\"\") in `METABOLITES` block (usually caused by excetrainious TAB at the end of line)."))

    if validate_features:
        if "MS_METABOLITE_DATA" in mwtabfile:
            from_metabolite_data_features = {i["metabolite_name"] for i in mwtabfile["MS_METABOLITE_DATA"]["MS_METABOLITE_DATA_START"]["DATA"]}
            if from_metabolites_features != from_metabolite_data_features:
                if "" in from_metabolite_data_features:
                    errors.append(ValueError(
                        "Feature with no name (\"\") in `MS_METABOLITE_DATA` block (usually caused by excetrainious "
                        "TAB at the end of line)."))
                if any(val for val in from_metabolite_data_features if val != ""):
                    errors.append(ValueError(
                        "`MS_METABOLITE_DATA` block contains additional features not found in `METABOLITES` block.\n\t "
                        "Additional features: {}".format(
                            from_metabolite_data_features.difference(from_metabolites_features))))

    return errors


def validate_file(mwtabfile, validate_samples=True, validate_factors=True, validate_features=True):
    """Validate entire ``mwTab`` formatted file one section at a time.

    :param mwtabfile: Instance of :class:`~mwtab.mwtab.MWTabFile`.
    :type mwtabfile: :class:`~mwtab.mwtab.MWTabFile`
    :param dict section_schema_mapping: Dictionary that provides mapping between section name and schema definition.
    :param validate_samples: Make sure that sample ids are consistent across file.
    :type validate_samples: :py:obj:`True` or :py:obj:`False`
    :param validate_factors: Make sure that factors are consistent across file.
    :type validate_factors: :py:obj:`True` or :py:obj:`False`
    :param validate_features: Make sure that features are consistent across file.
    :type validate_features: :py:obj:`True` or :py:obj:`False`
    :return: Validated file.
    :rtype: :py:class:`collections.OrderedDict`
    """
    errors = []
    try:
        errors.extend(_validate_samples_factors(mwtabfile, validate_samples, validate_factors))
        if mwtabfile.get("METABOLITES"):
            errors.extend(_validate_metabolites(mwtabfile, validate_features))
    except Exception:
        raise

    # for section_key, section in mwtabfile.items():
    #     try:
    #         schema = section_schema_mapping[section_key]
    #         section = validate_section(section=section, schema=schema)
    #         validated_mwtabfile[section_key] = section
    #     except Exception:
    #         raise

    return errors


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
    (r"(?i)(Human Metabolome D)[\S]{,}", "hmdb_id"),
]

duplicate_fields = {f: dict() for r, f in REGEXS}

if __name__ == '__main__':

    error_files = dict()
    unique_fields = dict()
    (_, _, filenames) = next(walk("/mlab/data/cdpo224/mwtab/data"))
    filenames = sorted(filenames)
    for filename in filenames:
        if not any(error in filename for error in processing_errors):
            mwfile = next(mwtab.read_files("/mlab/data/cdpo224/mwtab/data/{}".format(filename)))
            if mwfile.get("METABOLITES"):
                if mwfile["METABOLITES"]["METABOLITES_START"].get("Fields"):
                    from_metabolites_fields = set(mwfile["METABOLITES"]["METABOLITES_START"]["Fields"])
                    for field in from_metabolites_fields:
                        if field in unique_fields.keys():
                            unique_fields[field] += 1
                        else:
                            unique_fields[field] = 1
    del unique_fields["metabolite_name"]

    print(len(unique_fields))
    sent = 0
    items = list(unique_fields.items())
    for k, v in items:
        for r, f in REGEXS:
            if match(r, k):
                duplicate_fields[f].update({k: v})
                del unique_fields[k]
                break

    features = sorted(unique_fields.items(), key=lambda x: x[0].lower(), reverse=True)
    max_len = max(len(f[0]) for f in features)
    print("{{}:<{}} {}".format("FEATURE", str(max_len), "INSTANCES"))
    for f, n in features:
        print("{{}:<{}} {}".format(str(f), str(max_len), str(n)))

    print()
    for k in duplicate_fields.keys():
        print(k)
        for f in duplicate_fields[k].keys():
            print("\t{}".format(f))
