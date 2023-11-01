import hpotk
import pytest


from genophenocorr.preprocessing import PhenopacketPatientCreator
from genophenocorr.preprocessing import configure_patient_creator, load_phenopacket


class TestPhenopacketPatientCreator:

    @pytest.fixture
    def hpo(self) -> hpotk.MinimalOntology:
        return hpotk.load_minimal_ontology('testingDefaults/hp.json')

    @pytest.fixture
    def phenopacket_patient_creator(self, hpo: hpotk.MinimalOntology) -> PhenopacketPatientCreator:
        return configure_patient_creator(hpo)

    @pytest.mark.skip('Skipping online test')
    def test_load_phenopacket(self, phenopacket_patient_creator: PhenopacketPatientCreator):
        pp = load_phenopacket('../docs/data/simple_cohort/PMID_36446582_KBG12.json')
        patient = phenopacket_patient_creator.create_patient(pp)
        print(patient)

