import os

import pytest


@pytest.fixture(scope='session')
def fpath_preprocessing_data_dir() -> str:
    parent = os.path.dirname(__file__)
    return os.path.join(parent, 'data')
