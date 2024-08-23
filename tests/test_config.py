import os
import pathlib

import pytest

from gpsea.config import get_cache_dir_path, CACHE_ENV, DEFAULT_CACHE_PATH


class TestCacheDir:

    @pytest.mark.skipif(
        condition=CACHE_ENV in os.environ,
        reason="We cannot test default behavior under non-default condition (a variable is set)",
    )
    def test_get_cache_dir(self):
        cd = get_cache_dir_path()

        assert cd.name == DEFAULT_CACHE_PATH.name

    def test_create_using_environment_variable(
        self,
        tmp_path: pathlib.Path,
    ):
        target = tmp_path / ".ou_yeah"
        assert not target.exists()

        os.environ[CACHE_ENV] = str(target)

        cd = get_cache_dir_path()

        assert cd == target
