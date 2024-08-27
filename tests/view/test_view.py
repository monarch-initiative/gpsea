import hpotk
import pytest

from gpsea.model import Cohort
from gpsea.view import CohortViewable


class TestCohortViewable:

    @pytest.fixture
    def cohort_viewable(
            self,
            hpo: hpotk.MinimalOntology,
    ) -> CohortViewable:
        return CohortViewable(
            hpo=hpo,
        )

    def test_process(
            self,
            cohort_viewable: CohortViewable,
            toy_cohort: Cohort,
    ):
        toy_transcript_id = "NM_123.1"
        html = cohort_viewable.process(cohort=toy_cohort, transcript_id=toy_transcript_id)

        # A dummy test for now.
        assert len(html) != 0

