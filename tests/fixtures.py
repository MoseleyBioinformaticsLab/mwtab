
import pathlib
import time
import shutil

import pytest

def create_tmp_dir():
    path = pathlib.Path("tests/example_data/tmp/")
    path.mkdir(parents=True)

def delete_tmp_dir():
    path = pathlib.Path("tests/example_data/tmp/")
    if path.exists():
        shutil.rmtree(path)
        time_to_wait=10
        time_counter = 0
        while path.exists():
            time.sleep(1)
            time_counter += 1
            if time_counter > time_to_wait:
                raise FileExistsError(path + " was not deleted within " + str(time_to_wait) + " seconds, so it is assumed that it won't be and something went wrong.")


@pytest.fixture()
def init_tmp_dir():
    delete_tmp_dir()
    create_tmp_dir()
    yield
    delete_tmp_dir()

@pytest.fixture(autouse=True)
def init_tmp_dir_auto():
    delete_tmp_dir()
    create_tmp_dir()
    yield
    delete_tmp_dir()

@pytest.fixture()
def teardown_module():
    delete_tmp_dir()
    yield
    delete_tmp_dir()

@pytest.fixture(autouse=True)
def teardown_module_auto():
    delete_tmp_dir()
    yield
    delete_tmp_dir()

@pytest.fixture()
def init_tmp_dir_no_clean():
    create_tmp_dir()

@pytest.fixture()
def cleanup_tmp_dir():
    yield
    delete_tmp_dir()












