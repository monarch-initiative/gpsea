import os

import hpotk
import pytest

from gpsea.model import Cohort
from gpsea.view import CohortViewer


@pytest.mark.skip("Just for manual testing and debugging")
class TestCohortViewer:

    @pytest.fixture
    def cohort_viewer(
        self,
        hpo: hpotk.MinimalOntology,
    ) -> CohortViewer:
        return CohortViewer(
            hpo=hpo,
        )

    def test_process(
        self,
        cohort_viewer: CohortViewer,
        suox_cohort: Cohort,
        suox_mane_tx_id: str,
    ):
        report = cohort_viewer.process(
            cohort=suox_cohort,
            transcript_id=suox_mane_tx_id,
        )

        with open(os.path.join("dev", "SUOX.cohort.html"), "w") as fh:
            report.write(fh)
