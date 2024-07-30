import json

import hpotk
import pytest

from genophenocorr.io import GenophenocorrJSONEncoder, GenophenocorrJSONDecoder
from genophenocorr.model import Cohort
from genophenocorr.preprocessing import configure_caching_cohort_creator, load_phenopacket_folder


def test_round_trip(suox_cohort: Cohort):
    dumped = json.dumps(suox_cohort, cls=GenophenocorrJSONEncoder, indent=2)
    decoded = json.loads(dumped, cls=GenophenocorrJSONDecoder)

    assert suox_cohort == decoded


@pytest.mark.skip("Run manually to regenerate `suox_cohort`")
def test_regenerate_cohort(
    fpath_suox_cohort: str,
    hpo: hpotk.MinimalOntology,
):
    """
    The test for regenerating the `SUOX.json` file based on a cohort of phenopackets.
    The test needs path to a folder with phenopacket JSON files (empty `str` below).

    Note, the test may need to be run multiple times if the ENSEMBL API times out.
    """
    fpath_suox_pp_dir = "/path/to/SUOX/phenopackets"

    cohort_creator = configure_caching_cohort_creator(hpo, timeout=30)
    cohort = load_phenopacket_folder(
        fpath_suox_pp_dir, cohort_creator, validation_policy="strict"
    )
    with open(fpath_suox_cohort, "w") as fh:
        json.dump(cohort, fh, cls=GenophenocorrJSONEncoder, indent=2)
