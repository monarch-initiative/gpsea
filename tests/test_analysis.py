import hpotk
import pytest

import pandas as pd

from genophenocorr.analysis import configure_cohort_analysis
from genophenocorr.analysis.predicate import BooleanPredicate
from genophenocorr.data import get_toy_cohort
from genophenocorr.model import Cohort, VariantEffect


@pytest.mark.skip(reason='Disabled unless explicitly enabled')
class TestCommunistCohortAnalysis:

    @pytest.fixture
    def toy_cohort(self) -> Cohort:
        return get_toy_cohort()

    def test_compare_by_variant_effect(self, toy_cohort: Cohort, toy_hpo: hpotk.MinimalOntology):
        pd.set_option('expand_frame_repr', False)
        cohort_analysis = configure_cohort_analysis(toy_cohort, toy_hpo)
        results = cohort_analysis.compare_by_variant_effect(VariantEffect.MISSENSE_VARIANT, 'NM_1234.5')
        summary = results.summarize(toy_hpo, BooleanPredicate.YES)
        print(summary)
