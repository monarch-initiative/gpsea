import hpotk
import pytest

from gpsea.analysis.clf import HpoClassifier, DiseasePresenceClassifier
from gpsea.model import Cohort, Patient


def find_patient(pat_id: str, cohort: Cohort) -> Patient:
    for pat in cohort.all_patients:
        if pat.patient_id == pat_id:
            return pat
    raise ValueError(f'Could not find patient {pat_id}')


class TestHpoPredicate:

    @pytest.mark.parametrize('curie, patient_id, expected',
                             # Patient "HetSingleVar" has Phenotypes:
                             # Measured and present - 'HP:0001166;HP:0002266',  # Arachnodactyly;Focal clonic seizure
                             # Measured but excluded - 'HP:0001257',  # Spasticity
                             [
                                 # Test exact match
                                 ('HP:0001166',  # Arachnodactyly
                                  'HetSingleVar',
                                  "Yes"),
                                 # Test inferred annotations
                                 ('HP:0001250',  # Seizure
                                  'HetSingleVar',
                                  "Yes"),
                                 # Test excluded feature
                                 ('HP:0001257',  # Spasticity
                                  'HetSingleVar',
                                  "No"),
                             ])
    def test_phenotype_predicate__present_or_excluded(
            self,
            toy_cohort: Cohort,
            hpo: hpotk.MinimalOntology,
            curie: str,
            patient_id: str,
            expected: str,
    ):
        patient = find_patient(patient_id, toy_cohort)
        term_id = hpotk.TermId.from_curie(curie)
        predicate = HpoClassifier(hpo=hpo, query=term_id)
        actual = predicate.test(patient)

        assert actual is not None
        assert actual.phenotype == term_id
        assert actual.category.name == expected

    def test_phenotype_predicate__unknown(
            self,
            toy_cohort: Cohort,
            hpo: hpotk.MinimalOntology,
    ):
        # Not Measured and not Observed - 'HP:0006280',  # Chronic pancreatitis
        patient = find_patient('HetSingleVar', toy_cohort)
        term_id = hpotk.TermId.from_curie('HP:0006280')
        predicate = HpoClassifier(hpo=hpo, query=term_id)
        actual = predicate.test(patient)

        assert actual is None


class TestDiseasePresencePredicate:

    @pytest.mark.parametrize(
        'patient_id, patient_category',
        [
            ('HetSingleVar', "Yes"),
            ('HomoVar', "No"),
        ])
    def test_disease_predicate(
        self,
        patient_id: str,
        patient_category: str,
        toy_cohort: Cohort,

    ):
        patient = find_patient(patient_id, toy_cohort)
        disease_id = hpotk.TermId.from_curie("OMIM:148050")
        predicate = DiseasePresenceClassifier(disease_id)
        actual = predicate.test(patient)

        assert actual is not None
        assert actual.phenotype == disease_id
        assert actual.category.name == patient_category
