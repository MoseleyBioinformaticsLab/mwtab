# TODO: This file should not be included in the public release
import mwtab
from os import walk

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
    :return: None
    :rtype: :py:obj:`None`
    """
    errors = dict()

    from_subject_samples = {i["local_sample_id"] for i in mwtabfile["SUBJECT_SAMPLE_FACTORS"]["SUBJECT_SAMPLE_FACTORS"]}
    from_subject_factors = {i["factors"] for i in mwtabfile["SUBJECT_SAMPLE_FACTORS"]["SUBJECT_SAMPLE_FACTORS"]}

    if "" in from_subject_samples:
        errors.update({"0": ""})
    if "" in from_subject_factors:
        errors.update({"1": ""})

    if validate_samples:

        if "MS_METABOLITE_DATA" in mwtabfile:
            try:
                from_metabolite_data_samples = set(mwtabfile["MS_METABOLITE_DATA"]["MS_METABOLITE_DATA_START"]["Samples"])
                if not from_metabolite_data_samples.issubset(from_subject_samples):
                    if "" in from_metabolite_data_samples:
                        errors.update({"2": ""})
                    if any(val for val in from_metabolite_data_samples if val != ""):
                        errors.update({"3": from_metabolite_data_samples.difference(from_subject_samples)})
            except KeyError:
                errors.update({"4": ""})

        if "NMR_BINNED_DATA" in mwtabfile:
            try:
                from_nmr_binned_data_samples = set(mwtabfile["NMR_BINNED_DATA"]["NMR_BINNED_DATA_START"]["Fields"][1:])
                if not from_nmr_binned_data_samples.issubset(from_subject_samples):
                    if "" in from_nmr_binned_data_samples:
                        errors.update({"5": ""})
                    if any(val for val in from_nmr_binned_data_samples if val != ""):
                        errors.update({"6": from_nmr_binned_data_samples.difference(from_subject_samples)})
            except KeyError:
                errors.update({"7": ""})

    if validate_factors:

        if "MS_METABOLITE_DATA" in mwtabfile:
            try:
                from_metabolite_data_factors = set(mwtabfile["MS_METABOLITE_DATA"]["MS_METABOLITE_DATA_START"]["Factors"])
                if not from_metabolite_data_factors.issubset(from_subject_factors):
                    if "" in from_metabolite_data_factors:
                        errors.update({"8": ""})
                    if any(val for val in from_metabolite_data_factors if val != ""):
                        errors.update({"9": from_metabolite_data_factors.difference(from_subject_factors)})
            except KeyError:
                errors.update({"10": ""})

    return errors


def _validate_metabolites(mwtabfile, validate_features=True):
    """Validate section of ``mwTab`` formatted file.

    :param section: Section of :class:`~mwtab.mwtab.MWTabFile`.
    :param schema: Schema definition.
    :return: Validated section.
    :rtype: :py:class:`collections.OrderedDict`
    """
    errors = dict()

    from_metabolites_features = {i["metabolite_name"] for i in mwtabfile["METABOLITES"]["METABOLITES_START"]["DATA"]}

    if validate_features:
        if "MS_METABOLITE_DATA" in mwtabfile:
            from_metabolite_data_features = {i["metabolite_name"] for i in mwtabfile["MS_METABOLITE_DATA"]["MS_METABOLITE_DATA_START"]["DATA"]}
            if from_metabolites_features != from_metabolite_data_features:
                print(mwtabfile.analysis_id)
                errors.update({"11": ""})

    schema_keys = {"metabolite_name", "moverz_quant", "ri", "ri_type", "pubchem_id", "inchi_key", "kegg_id", "other_id", "other_id_type"}
    if mwtabfile["METABOLITES"]["METABOLITES_START"].get("Fields"):
        section_keys = set(mwtabfile["METABOLITES"]["METABOLITES_START"]["Fields"])
    else:
        errors.update({"15": ""})
        section_keys = set(mwtabfile["METABOLITES"]["METABOLITES_START"]["DATA"][0].keys())

    diff = section_keys.difference(schema_keys)
    if diff:
        return {'14': diff}


def validate_file(mwtabfile, validate_samples=True, validate_factors=True):
    """Validate entire ``mwTab`` formatted file one section at a time.

    :param mwtabfile: Instance of :class:`~mwtab.mwtab.MWTabFile`.
    :type mwtabfile: :class:`~mwtab.mwtab.MWTabFile`
    :param dict section_schema_mapping: Dictionary that provides mapping between section name and schema definition.
    :param validate_samples: Make sure that sample ids are consistent across file.
    :type validate_samples: :py:obj:`True` or :py:obj:`False`
    :param validate_factors: Make sure that factors are consistent across file.
    :type validate_factors: :py:obj:`True` or :py:obj:`False`
    :return: Validated file.
    :rtype: :py:class:`collections.OrderedDict`
    """
    errors = dict()
    errors.update(_validate_samples_factors(mwtabfile, validate_samples, validate_factors))

    if mwtabfile.get("METABOLITES"):
        print(mwtabfile.analysis_id)
        metabolites_error = _validate_metabolites(mwtabfile)
        if metabolites_error:
            errors.update(metabolites_error)

    return errors


if __name__ == '__main__':

    error_files = dict()
    (_, _, filenames) = next(walk("/mlab/data/cdpo224/mwtab/data"))
    filenames = sorted(filenames)
    for filename in filenames:
        if not any(error in filename for error in processing_errors):
            mwfile = next(mwtab.read_files("/mlab/data/cdpo224/mwtab/data/{}".format(filename)))
            errors = validate_file(mwfile)
            if errors:
                error_files.update({filename[:-4]: errors})

    fields = set()
    for fn in error_files.keys():
        if error_files[fn].get("14"):
            fields.update(error_files[fn]["14"])
    print(fields)
    print(len(fields))
