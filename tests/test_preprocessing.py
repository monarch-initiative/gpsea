import hpotk
import pytest


from genophenocorr.preprocessing import PhenopacketPatientCreator, PhenotypeCreator
from genophenocorr.preprocessing import configure_patient_creator, load_phenopacket


class TestPhenopacketPatientCreator:

    @pytest.fixture
    def phenopacket_patient_creator(self, toy_hpo: hpotk.MinimalOntology) -> PhenopacketPatientCreator:
        return configure_patient_creator(toy_hpo)

    @pytest.mark.skip('Skipping online test')
    def test_load_phenopacket(self, phenopacket_patient_creator: PhenopacketPatientCreator):
        pp = load_phenopacket('../docs/data/simple_cohort/PMID_36446582_KBG12.json')
        patient = phenopacket_patient_creator.create_patient(pp)
        print(patient)


class TestPhenotypeCreator:

    @pytest.fixture
    def creator(self, toy_hpo: hpotk.MinimalOntology,
                toy_validation_runner: hpotk.validate.ValidationRunner) -> PhenotypeCreator:
        return PhenotypeCreator(toy_hpo, toy_validation_runner)

    def test_something(self, creator: PhenotypeCreator):
        inputs = (
            ('HP:0001250', True),
        )
        out = creator.process(inputs)
        print(out)
