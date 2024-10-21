import io

import hpotk
import pytest

from gpsea.model import Cohort
from gpsea.view import CohortViewer


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
        toy_cohort: Cohort,
    ):
        toy_transcript_id = "NM_123.1"

        report = cohort_viewer.process(
            cohort=toy_cohort, transcript_id=toy_transcript_id
        )

        buf = io.StringIO()
        report.write(buf)
        html = buf.getvalue()

        # A dummy test for now.
        assert len(html) != 0
