from os import walk
import re
import mwtab
from mwtab import read_files, validate_file
# section_schema_mapping

processing_errors = [
    "AN000404", "AN000405", "AN000410", "AN000415", "AN000436", "AN000439", "AN001311", "AN001312", "AN001313",
    "AN000696", "AN001467", "AN001576", "AN001684", "AN001685", "AN001721", "AN001979", "AN002035", "AN001492",
    "AN001493", "AN001499", "AN001761", "AN001762", "AN001776", "AN001777", "AN001982", "AN001689", "AN001690",
    "AN001992",

    "AN000144", "AN000168", "AN000172", "AN000416", "AN000603", "AN000611", "AN000612", "AN000615", "AN000616",
    "AN000617", "AN000632", "AN001069"

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
    (r"(?i)(Human Metabolome D)[\S]{,}", "hmdb_id"),
]

duplicate_fields = {f: dict() for r, f in REGEXS}

if __name__ == '__main__':

    count = 0
    (_, _, filenames) = next(walk("/mlab/data/cdpo224/mwtab/data"))
    filenames = sorted(filenames)
    for filename in filenames:
        try:
            mwfile = next(read_files("/mlab/data/cdpo224/mwtab/data/{}".format(filename)))
            if re.search(r"(?i)(human)", mwfile["SUBJECT"]["SUBJECT_TYPE"]):
                count += 1
        except Exception:
            pass

    print(count)
    exit()


    error_files = {
        "processing": set(),
        "validating": set()
    }
    unique_fields = dict()
    (_, _, filenames) = next(walk("/mlab/data/cdpo224/mwtab/data"))
    filenames = sorted(filenames)
    for filename in filenames:
        print(filename.split(".")[0])
        try:
            mwfile = next(read_files("/mlab/data/cdpo224/mwtab/data/{}".format(filename)))

            try:
                errors = validate_file(mwfile, test=True)
                if errors:
                    error_files["validating"].add(filename.split(".")[0])
                    # print("\tValidating Error(s)")
                    # for e in errors:
                    #     print("\t\t" + e.__class__.__name__)
                    #     print("\t\t" + str(e))
            except KeyError as e:
                error_files["validating"].add(filename.split(".")[0])
                print("\tValidating Error(s)")
                print("\t\t" + e.__class__.__name__)
                print("\t\t" + str(e))

        except Exception as e:
            error_files["processing"].add(filename.split(".")[0])
            # print("\tProcessing Error")
            # print("\t\t"+e.__class__.__name__)
            # print("\t\t"+str(e))

    print(len(error_files["validating"]))
    print(len(error_files["processing"]))

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
