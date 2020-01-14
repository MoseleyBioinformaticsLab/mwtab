import mwtab
import re
import mwtab.validator as mwval
import mwtab.mwrest as mwr
import jsonpickle
from urllib.request import urlopen
from os import walk
TESTYPATH = '/mlab/data/cdpo224/mwtab/testymctestyface/{}.txt'
DATAPATH = '/mlab/data/cdpo224/mwtab/data/'

mwtab.fileio.VERBOSE = True


# Useful ReGeX r'(AN[0-9]{6})'
processing_errors = [
    "AN000404", "AN000405", "AN000410", "AN000415", "AN000436", "AN000439", "AN001311", "AN001312", "AN001313",
    "AN000696", "AN001467", "AN001576", "AN001684", "AN001685", "AN001721", "AN001979", "AN002035", "AN001492",
    "AN001493", "AN001499", "AN001761", "AN001762", "AN001776", "AN001777", "AN001982", "AN001689", "AN001690",
    "AN001992"
]


def gnerate_filenames(data_path="data/"):
    (_, _, filenames) = next(walk(DATAPATH))
    for filename in filenames:
        yield filename


def access_single_file(analysis_id):
    """
    Access Single File
    """
    mwfile = next(mwtab.read_files(analysis_id))
    return mwfile


def access_saved(filenames):
    """
    Access Each File in `data/`
    """
    for filename in filenames:
        with open(DATAPATH + str(filename), 'r') as f:
            yield f
            f.close()


def open_error_file():
    with open("errors.json", 'r') as error_file:
        json_str = error_file.read()
        error_file.close()

    return jsonpickle.decode(json_str)


def fetch_all():
    """
    Pull Down Every Analysis
    """
    analysis_ids = mwr.analysis_ids()
    urls = mwr.generate_mwtab_urls(*analysis_ids)

    for url in urls:
        print(url)
        response = urlopen(url)
        path = response.read().decode('utf-8')
        response.close()
        with open("data/{}.txt".format(re.search(r'(AN[0-9]{6})', url).group(0)), "w") as outfile:
            outfile.write(path)
            outfile.close()


if __name__ == '__main__':

    mwfile = next(mwtab.read_files("data/AN000009.txt"))
    from_metabolites_features = {i["metabolite_name"] for i in mwfile["METABOLITES"]["METABOLITES_START"]["DATA"]}
    from_metabolite_data_features = {feature["metabolite_name"] for feature in mwfile["MS_METABOLITE_DATA"]["MS_METABOLITE_DATA_START"]["DATA"]}
    print(from_metabolite_data_features)
    exit()

    """
    Check every analysis for errors
    """
    errors = list()
    (_, _, filenames) = next(walk("/mlab/data/cdpo224/mwtab/data"))
    filenames = sorted(filenames)
    for filename in filenames:
        if not any(error in filename for error in processing_errors):
            print(filename)
            try:
                mwfile = next(mwtab.read_files("data/{}".format(filename), validate=True))
            except Exception as e:
                print(e)
                errors.append(re.search(r'(AN[0-9]{6})', filename).group(0))
            print()

    print(len(errors))
    with open("errors.json", "w") as outfile:
        outfile.write(jsonpickle.encode(errors))
        outfile.close()
