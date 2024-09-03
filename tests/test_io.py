import json
import pathlib

import hpotk
import pytest

from gpsea.io import GpseaJSONEncoder, GpseaJSONDecoder
from gpsea.model import Cohort
from gpsea.preprocessing import configure_caching_cohort_creator, load_phenopackets


def test_round_trip(suox_cohort: Cohort):
    dumped = json.dumps(suox_cohort, cls=GpseaJSONEncoder, indent=2)
    decoded = json.loads(dumped, cls=GpseaJSONDecoder)

    assert suox_cohort == decoded


@pytest.mark.skip("Run manually to regenerate `suox_cohort`")
def test_regenerate_cohort(
    fpath_suox_cohort: str,
    hpo: hpotk.MinimalOntology,
    tmp_path: pathlib.Path,
):
    """
    The test for regenerating the `SUOX.json` file based on a cohort of phenopackets.

    Note, the test may need to be run multiple times if the ENSEMBL API times out.
    """
    from ppktstore.registry import configure_phenopacket_registry

    registry = configure_phenopacket_registry(store_dir=tmp_path)
    with registry.open_phenopacket_store('0.1.18') as ps:
        phenopackets = tuple(ps.iter_cohort_phenopackets('SUOX'))
    
    cohort_creator = configure_caching_cohort_creator(hpo, timeout=30.)
    cohort, validation = load_phenopackets(
        phenopackets=phenopackets,
        cohort_creator=cohort_creator,
        validation_policy="strict",
    )

    if not validation.is_ok():
        raise ValueError('The cohort MUST be OK!')

    with open(fpath_suox_cohort, "w") as fh:
        json.dump(cohort, fh, cls=GpseaJSONEncoder, indent=2)
