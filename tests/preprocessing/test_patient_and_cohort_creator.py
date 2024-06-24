import hpotk
import pytest


from genophenocorr.preprocessing import PhenopacketPatientCreator, CohortCreator
from genophenocorr.preprocessing import configure_patient_creator, configure_cohort_creator, load_phenopacket
from genophenocorr.preprocessing import load_phenopacket_folder


class TestPhenopacketPatientCreator:

    @pytest.fixture
    def phenopacket_patient_creator(self, hpo: hpotk.MinimalOntology) -> PhenopacketPatientCreator:
        return configure_patient_creator(hpo)

    @pytest.mark.skip('Skipping online test')
    def test_load_phenopacket(self, phenopacket_patient_creator: PhenopacketPatientCreator):
        pp = load_phenopacket('../docs/data/simple_cohort/PMID_36446582_KBG12.json')
        patient = phenopacket_patient_creator.create_patient(pp)
        print(patient)


class TestPhenopacketCohortCreator:

    @pytest.fixture
    def phenopacket_cohort_creator(self, hpo: hpotk.MinimalOntology) -> CohortCreator:
        return configure_cohort_creator(hpo)

    @pytest.mark.skip('Skipping online test')
    def test_load_phenopacket(self, phenopacket_cohort_creator: CohortCreator):
        cohort = load_phenopacket_folder('docs/data/simple_cohort', phenopacket_cohort_creator)
        print(cohort)
