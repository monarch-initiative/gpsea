import os
import pathlib
import typing


CACHE_ENV = "GPSEA_CACHEDIR"
"""
The name of the environment variable consulted by GPSEA
to set the cache directory.
"""

DEFAULT_CACHE_PATH = pathlib.Path(".gpsea_cache")
"""
Default path to GPSEA cache directory.
"""


def get_cache_dir_path(
    cache: typing.Optional[typing.Union[str, pathlib.Path]] = None,
) -> pathlib.Path:
    """
    Get path to the default cache directory.
    
    First try to use `cache` argument.
    If `cache` is `None`, then use the `GPSEA_CACHE` environment variable, if set.
    Last, fall back to default cache path (`.gpsea_cache` in the current working directory).
    
    Note: the cache directory is *not* created if it does not exist.

    :path cache: a `str` or a :class:`~pathlib.Path` to fallback cache
        or `None` if the default cache folder should be used.
    """
    if cache is None:
        if CACHE_ENV in os.environ:
            # Let's use the environment variable.
            cache = os.environ[CACHE_ENV]
        else:
            # Nothing provided as environment variable, let's use the default.
            cache = DEFAULT_CACHE_PATH
    
    # Ensure `cache_path` is path
    if isinstance(cache, pathlib.Path):
        cache_path = cache
    elif isinstance(cache, str):
        cache_path = pathlib.Path(cache)
    else:
        raise ValueError(
            f"`cache` must be a `str` or `pathlib.Path` but was {type(cache)}"
        )

    return cache_path
