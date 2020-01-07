# TODO: This file should not be included in the public release
import mwtab
from os import walk

processing_errors = [
    "AN000404", "AN000405", "AN000410", "AN000415", "AN000436", "AN000439", "AN001311", "AN001312", "AN001313",
    "AN000696", "AN001467", "AN001576", "AN001684", "AN001685", "AN001721", "AN001979", "AN002035", "AN001492",
    "AN001493", "AN001499", "AN001761", "AN001762", "AN001776", "AN001777", "AN001982", "AN001689", "AN001690",
    "AN001992"
]

if __name__ == '__main__':

    has_mnr_data = list()
    has_ms_data = list()
    has_raw_file = list()

    (_, _, filenames) = next(walk("/mlab/data/cdpo224/mwtab/data"))
    filenames = sorted(filenames)
    for filename in filenames:
        if not any(error in filename for error in processing_errors):
            mwfile = next(mwtab.read_files("/mlab/data/cdpo224/mwtab/data/{}".format(filename)))
