import os

import hpotk
import pytest

from gpsea.analysis.pcats import HpoTermAnalysisResult
from gpsea.model import Cohort
from gpsea.view import (
    CohortViewer,
    CohortVariantViewer,
    MtcStatsViewer,
)


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


@pytest.mark.skip("Just for manual testing and debugging")
def test_viewer(
    suox_mane_tx_id: str,
    suox_cohort: Cohort,
):
    viewer = CohortVariantViewer(tx_id=suox_mane_tx_id)
    html = viewer.process(suox_cohort)

    with open("all_variants.html", "w") as fh:
        html.write(fh)


class TestMtcStatsViewer:

    @pytest.fixture
    def stats_viewer(self) -> MtcStatsViewer:
        return MtcStatsViewer()

    @pytest.mark.skip("Just for manual testing and debugging")
    def test_process(
        self,
        stats_viewer: MtcStatsViewer,
        hpo_term_analysis_result: HpoTermAnalysisResult,
    ):
        report = stats_viewer.process(result=hpo_term_analysis_result)
        with open("mtc_stats.html", "w") as fh:
            report.write(fh)
