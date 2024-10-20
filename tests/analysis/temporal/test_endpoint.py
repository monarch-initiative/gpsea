import hpotk
import pytest

from gpsea.model import Patient, Age, VitalStatus, Status, Disease
from gpsea.analysis.temporal.endpoint import death, disease_onset, hpo_onset
from gpsea.model import Phenotype


@pytest.fixture
def patient_no_data() -> Patient:
    return Patient.from_raw_parts(
        labels="?",
    )


@pytest.fixture
def alive() -> Patient:
    return Patient.from_raw_parts(
        labels="A",
        age=Age.postnatal_days(days=40),
        vital_status=VitalStatus(status=Status.ALIVE, age_of_death=None),
    )


@pytest.fixture
def deceased() -> Patient:
    return Patient.from_raw_parts(
        labels="D",
        age=Age.postnatal_days(days=60),
        vital_status=VitalStatus(
            status=Status.DECEASED, age_of_death=Age.postnatal_days(60)
        ),
    )


class TestDeath:

    def test_compute_survival__alive(
        self,
        alive: Patient,
    ):
        endpoint = death(kind="postnatal")

        survival = endpoint.compute_survival(alive)

        assert survival is not None
        assert survival.is_censored
        assert survival.value == pytest.approx(40.0)

    def test_compute_survival__deceased(
        self,
        deceased: Patient,
    ):
        endpoint = death(kind="postnatal")

        survival = endpoint.compute_survival(deceased)

        assert survival is not None
        assert not survival.is_censored
        assert survival.value == pytest.approx(60.0)

    def test_compute_survival__unknown(
        self,
        patient_no_data: Patient,
    ):
        endpoint = death(kind="postnatal")

        survival = endpoint.compute_survival(patient_no_data)

        assert survival is None

    def test_display_question(self):
        endpoint = death(kind="postnatal")

        assert endpoint.display_question() == "Compute time until postnatal death"


class TestDiseaseOnset:
    
    @pytest.fixture(scope="class")
    def with_disease(self):
        return Patient.from_raw_parts(
            labels="X",
            diseases=(
                Disease.from_raw_parts(
                    term_id="OMIM:123456",
                    name="NOT REALLY IMPORTANT HERE",
                    is_observed=True,
                    onset=Age.postnatal_days(40),
                ),
            ),
        )

    def test_compute_survival__alive(
        self,
        alive: Patient,
    ):
        endpoint = disease_onset(disease_id="OMIM:123456")

        survival = endpoint.compute_survival(alive)

        assert survival is not None
        assert survival.is_censored
        assert survival.value == pytest.approx(40.0)

    def test_compute_survival__with_target_disease(
        self,
        with_disease: Patient,
    ):
        endpoint = disease_onset(disease_id="OMIM:123456")
        
        survival = endpoint.compute_survival(with_disease)
        
        assert survival is not None
        assert not survival.is_censored
        assert survival.value == pytest.approx(40.)
    
    def test_compute_survival__off_target_disease(
        self,
        with_disease: Patient,
    ):
        endpoint = disease_onset(disease_id="OMIM:100000")

        survival = endpoint.compute_survival(with_disease)

        assert survival is None

    def test_display_question(self):
        endpoint = disease_onset(disease_id="OMIM:100000")

        assert endpoint.display_question() == "Compute time until postnatal diagnosis of OMIM:100000"


class TestPhenotypeOnset:

    @pytest.fixture
    def phenotyped(self) -> Patient:
        return Patient.from_raw_parts(
            labels="P",
            age=Age.postnatal_years(30),
            vital_status=VitalStatus(status=Status.ALIVE, age_of_death=None),
            phenotypes=(
                Phenotype.from_raw_parts(
                    "HP:0007359",  # Focal-onset seizure
                    is_observed=True,
                    onset=Age.postnatal_days(20),
                ),
                Phenotype.from_raw_parts(
                    "HP:0020219",  # Motor seizure
                    is_observed=True,
                    onset=Age.postnatal_days(25),
                ),
                Phenotype.from_raw_parts(
                    "HP:0011153",  # Focal motor seizure
                    is_observed=False,
                    onset=Age.postnatal_days(25),
                )
            )
        )

    def test_compute_survival__alive(
        self,
        hpo: hpotk.MinimalOntology,
        alive: Patient,
    ):
        endpoint = hpo_onset(hpo, term_id="HP:0001250")  # Seizure

        survival = endpoint.compute_survival(alive)
        
        assert survival is not None
        assert survival.is_censored
        assert survival.value == pytest.approx(40.)

    def test_compute_survival__phenotyped(
        self,
        hpo: hpotk.MinimalOntology,
        phenotyped: Patient,
    ):
        endpoint = hpo_onset(hpo, term_id="HP:0001250")  # Seizure

        survival = endpoint.compute_survival(phenotyped)

        assert survival is not None
        assert not survival.is_censored
        assert survival.value == pytest.approx(20.)  # Focal-onset seizure has the earliest onset

    def test_display_question(
        self,
        hpo: hpotk.MinimalOntology,
    ):
        endpoint = hpo_onset(hpo, term_id="HP:0001250")  # Seizure

        assert endpoint.display_question() == "Compute time until postnatal onset of Seizure"
