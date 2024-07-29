import hpotk
import pytest


from genophenocorr.preprocessing import configure_caching_cohort_creator, CohortCreator, load_phenopacket_folder


class TestPhenopacketCohortCreator:

    @pytest.fixture
    def phenopacket_cohort_creator(self, hpo: hpotk.MinimalOntology) -> CohortCreator:
        return configure_caching_cohort_creator(hpo)

    @pytest.mark.skip('Skipping online test')
    def test_load_phenopacket(self, phenopacket_cohort_creator: CohortCreator):
        cohort = load_phenopacket_folder('docs/data/simple_cohort', phenopacket_cohort_creator)
        print(cohort)
