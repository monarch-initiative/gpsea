import pytest

from gpsea.model import Cohort
from gpsea.view import CohortVariantViewer


@pytest.mark.skip("For manual run only")
def test_viewer(
    suox_mane_tx_id: str,
    suox_cohort: Cohort,
):
    viewer = CohortVariantViewer(transcript_id=suox_mane_tx_id)
    html = viewer.process(suox_cohort)

    with open("all_variants.html", "w") as fh:
        fh.write(html)
