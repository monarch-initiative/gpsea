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

    def test_process_suox_cohort(
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

    def test_process_cyp21a2_cohort(
        self,
        cohort_viewer: CohortViewer,
        cyp21a2_cohort: Cohort,
        cyp21a2_mane_tx_id: str,
    ):
        report = cohort_viewer.process(
            cohort=cyp21a2_cohort,
            transcript_id=cyp21a2_mane_tx_id,
        )

        with open(os.path.join("dev", "CYP21A2.cohort.html"), "w") as fh:
            report.write(fh)
