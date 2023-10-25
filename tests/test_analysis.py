import hpotk
import pytest

from genophenocorr.analysis import CommunistCohortAnalysis
from genophenocorr.data import get_toy_cohort
from genophenocorr.model import Cohort, VariantEffect


# @pytest.mark.skip(reason='Disabled unUnless explicitly enabled')
class TestSlimCohortAnalysis:

    @pytest.fixture
    def toy_cohort(self) -> Cohort:
        return get_toy_cohort()

    @pytest.fixture
    def hpo(self) -> hpotk.MinimalOntology:
        return hpotk.load_minimal_ontology('')

    def test_compare_by_variant_effect(self, toy_cohort: Cohort, hpo: hpotk.MinimalOntology):
        cohort_analysis = CommunistCohortAnalysis(toy_cohort, hpo)
        results = cohort_analysis.compare_by_variant_effect(VariantEffect.MISSENSE_VARIANT, 'NM_1234.5')
        print(results)
