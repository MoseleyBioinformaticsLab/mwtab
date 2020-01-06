import mwtab
from os import walk

"""
Error List
0 - Missing Sample ID in `SUBJECT_SAMPLE_FACTORS` block. (should be impossible...)
1 - Extra Sample ID in `SUBJECT_SAMPLE_FACTORS` block.
2 - Extra tab in `MS_METABOLITE_DATA` block (samples).
3 - Extra Sample ID in `MS_METABOLITE_DATA` block.

4 - Extra tab in `NMR_BINNED_DATA` block.
5 - Extra Sample ID in `NMR_BINNED_DATA` block.

6 - Missing Sample Factor in `SUBJECT_SAMPLE_FACTORS` block. (should be impossible...)
7 - Extra Sample Factor in `SUBJECT_SAMPLE_FACTORS` block.
8 - Extra tab in `MS_METABOLITE_DATA` block (factors).
9 - Extra Factor in `MS_METABOLITE_DATA` block.
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
    from_subject_samples = {i["local_sample_id"] for i in mwtabfile["SUBJECT_SAMPLE_FACTORS"]["SUBJECT_SAMPLE_FACTORS"]}
    from_subject_factors = {i["factors"] for i in mwtabfile["SUBJECT_SAMPLE_FACTORS"]["SUBJECT_SAMPLE_FACTORS"]}

    if validate_samples:

        if "MS_METABOLITE_DATA" in mwtabfile:
            from_metabolite_data_samples = set(mwtabfile["MS_METABOLITE_DATA"]["MS_METABOLITE_DATA_START"]["Samples"])
            if from_subject_samples != from_metabolite_data_samples:
                from_ss = from_subject_samples.difference(from_metabolite_data_samples)
                from_mds = from_metabolite_data_samples.difference(from_subject_samples)
                if from_ss:
                    if from_ss == {''}:
                        return {'0': ''}
                    elif '' in from_ss:
                        return {'0': '', '1': from_ss}
                    else:
                        return {'1': from_ss}
                if from_mds:
                    if from_mds == {''}:
                        return {'2': ''}
                    elif '' in from_mds:
                        return {'2': '', '3': from_mds}
                    else:
                        return {'3': from_mds}

        if "NMR_BINNED_DATA" in mwtabfile:
            from_nmr_binned_data_samples = set(mwtabfile["NMR_BINNED_DATA"]["NMR_BINNED_DATA_START"]["Fields"][1:])
            if from_subject_samples != from_nmr_binned_data_samples:
                from_ss = from_subject_samples.difference(from_nmr_binned_data_samples)
                from_nbds = from_nmr_binned_data_samples.difference(from_subject_samples)
                if from_ss:
                    if from_ss == {''}:
                        return {'0': ''}
                    elif '' in from_ss:
                        return {'0': '', '1': from_ss}
                    else:
                        return {'1': from_ss}
                if from_nbds:
                    if from_nbds == {''}:
                        return {'4': ''}
                    elif '' in from_nbds:
                        return {'4': '', '5': from_nbds}
                    else:
                        return {'5': from_nbds}

    if validate_factors:

        if "MS_METABOLITE_DATA" in mwtabfile:
            from_metabolite_data_factors = set(mwtabfile["MS_METABOLITE_DATA"]["MS_METABOLITE_DATA_START"]["Factors"])
            if from_subject_factors != from_metabolite_data_factors:
                from_sf = from_subject_factors.difference(from_metabolite_data_factors)
                from_mdf = from_metabolite_data_factors.difference(from_subject_factors)
                if from_sf:
                    if from_sf == {''}:
                        return {'6': ''}
                    elif '' in from_sf:
                        return {'6': '', '7': from_sf}
                    else:
                        return {'7': from_sf}
                if from_mdf:
                    if from_mdf == {''}:
                        return {'8': ''}
                    elif '' in from_mdf:
                        return {'8': '', '9': from_mdf}
                    else:
                        return {'9': from_mdf}


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
    errors = _validate_samples_factors(mwtabfile, validate_samples, validate_factors)

    return errors


if __name__ == '__main__':
    (_, _, filenames) = next(walk("/mlab/data/cdpo224/mwtab/data"))
    filenames = sorted(filenames)
    for filename in filenames:
        if not any(error in filename for error in processing_errors):
            print(filename)
            mwfile = next(mwtab.read_files("/mlab/data/cdpo224/mwtab/data/{}".format(filename)))
            print(validate_file(mwfile))
