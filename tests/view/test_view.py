import hpotk
import pytest

from genophenocorr.model import Cohort
from genophenocorr.view import CohortViewable


class TestCohortViewable:

    @pytest.fixture
    def cohort_viewable(
            self,
            toy_hpo: hpotk.MinimalOntology,
    ) -> CohortViewable:
        return CohortViewable(
            hpo=toy_hpo,
        )

    def test_process(
            self,
            cohort_viewable: CohortViewable,
            toy_cohort: Cohort,
    ):
        html = cohort_viewable.process(cohort=toy_cohort)

        # A dummy test for now.
        assert len(html) != 0

