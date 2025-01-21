import typing

from gpsea.preprocessing import configure_default_protein_metadata_service

from ._base import BaseProteinVisualizer, CohortArtist
from ._protein_visualizer import ProteinVisualizer

def configure_default_protein_visualizer(
    random_seed: int = 42,
) -> BaseProteinVisualizer:
    return ProteinVisualizer(
        random_seed=random_seed,
    )

def configure_default_cohort_artist(
    protein_source: typing.Literal["UNIPROT"] = "UNIPROT",
    cache_dir: typing.Optional[str] = None,
    timeout: typing.Union[float, int] = 30.0,
    random_seed: int = 42,
) -> CohortArtist:
    """
    Get the default :class:`~gpsea.view.CohortArtist`.
    
    :param protein_source: a `str` with the code of the protein data sources (currently accepting just `UNIPROT`).
    :param cache_dir: path to the folder where we will cache the results fetched from the remote APIs or `None`
        if the data should be cached as described by :func:`~gpsea.config.get_cache_dir_path` function.
        In any case, the directory will be created if it does not exist (including any non-existing parents).
    :param timeout: a `float` or an `int` for the timeout in seconds for the REST APIs.
    """
    protein_visualizer = configure_default_protein_visualizer(
        random_seed=random_seed,
    )
    protein_metadata_service = configure_default_protein_metadata_service(
        protein_source=protein_source,
        cache_dir=cache_dir,
        timeout=timeout,
    )
    return CohortArtist(
        protein_visualizer=protein_visualizer,
        protein_metadata_service=protein_metadata_service,
    )
