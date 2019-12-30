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

    #         f.close()
    #     if '_RESULTS_FILE' in mwtab_str:
    #         include_results_file.append(filename.split('.')[0])
    #
    # print(sorted(include_results_file))
    # print(len(include_results_file))


def open_error_file():
    with open("errors.json", 'r') as error_file:
        json_str = error_file.read()
        error_file.close()

    return jsonpickle.decode(json_str)


if __name__ == '__main__':

    # 1878

    access_single_file('917')
    exit()
    """
    Pull Down Every Analysis
    """
    # analysis_ids = mwr.analysis_ids()
    # urls = mwr.generate_mwtab_urls(*analysis_ids)
    #
    # for url in urls:
    #     print(url)
    #     response = urlopen(url)
    #     path = response.read().decode('utf-8')
    #     response.close()
    #     with open("data/{}.txt".format(re.search(r'(AN[0-9]{6})', url).group(0)), "w") as outfile:
    #         outfile.write(path)

    # for mwfile in mwtab.read_files("data/AN001818.txt"):
    #     print(mwfile.analysis_id)
    # exit()


    #
    # errors["Processing Error"]["Format"] = sorted(list(errors["Processing Error"]["Format"]))
    #
    # for an_id in errors["Processing Error"]["Format"]:
    #     print(an_id)
    #     try:
    #         mwfile = next(mwtab.read_files("data/{}.txt".format(an_id)))
    #     except Exception:
    #         pass

    # error_ids = set(errors["Processing Error"])
    #
    # errors = {"Processing Error": {"Blank": set(), "Format": set()}}
    # blanks = list()
    #
    # for error_id in error_ids:
    #     with open('data/{}.txt'.format(error_id), 'r') as f:
    #         file_str = f.read()
    #         f.close()
    #     if file_str == '':
    #         errors["Processing Error"]["Blank"].add(error_id)
    #
    # errors["Processing Error"]["Format"] = error_ids.difference(errors["Processing Error"]["Blank"])
    #
    # with open("errors.json", 'w') as error_file:
    #     error_file.write(jsonpickle.encode(errors))
    #     error_file.close()

    """
    Check every analysis for errors
    """
    errors = list()
    (_, _, filenames) = next(walk("/mlab/data/cdpo224/mwtab/data"))
    filenames = sorted(filenames)
    for filename in filenames:
        print(filename)
        try:
            mwfile = next(mwtab.read_files("data/{}".format(filename)))
        except Exception as e:
            print(e)
            errors.append(re.search(r'(AN[0-9]{6})', filename).group(0))
        print()

    errors = {'Processing Error': errors}
