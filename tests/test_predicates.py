import hpotk
import pytest

from gpsea.analysis.predicate import PatientCategory, PatientCategories
from gpsea.analysis.predicate.phenotype import PropagatingPhenotypePredicate, DiseasePresencePredicate
from gpsea.analysis.predicate.genotype import *
from gpsea.model import Cohort, Patient, FeatureType, VariantEffect
from gpsea.model.genome import Region
from gpsea.preprocessing import ProteinMetadataService


def find_patient(pat_id: str, cohort: Cohort) -> Patient:
    for pat in cohort.all_patients:
        if pat.patient_id == pat_id:
            return pat
    raise ValueError(f'Could not find patient {pat_id}')


class TestPropagatingPhenotypeBooleanPredicate:

    @pytest.mark.parametrize('curie, patient_id, expected',
                             # Patient "HetSingleVar" has Phenotypes:
                             # Measured and present - 'HP:0001166;HP:0002266',  # Arachnodactyly;Focal clonic seizure
                             # Measured but excluded - 'HP:0001257',  # Spasticity
                             [
                                 # Test exact match
                                 ('HP:0001166',  # Arachnodactyly
                                  'HetSingleVar',
                                  PatientCategories.YES),
                                 # Test inferred annotations
                                 ('HP:0001250',  # Seizure
                                  'HetSingleVar',
                                  PatientCategories.YES),
                                 # Test excluded feature
                                 ('HP:0001257',  # Spasticity
                                  'HetSingleVar',
                                  PatientCategories.NO),
                             ])
    def test_phenotype_predicate__present_or_excluded(
            self,
            toy_cohort: Cohort,
            hpo: hpotk.MinimalOntology,
            curie: str,
            patient_id: str,
            expected: PatientCategory,
    ):
        patient = find_patient(patient_id, toy_cohort)
        term_id = hpotk.TermId.from_curie(curie)
        predicate = PropagatingPhenotypePredicate(hpo=hpo, query=term_id)
        actual = predicate.test(patient)

        assert actual.phenotype == term_id
        assert actual.category == expected

    def test_phenotype_predicate__unknown(
            self,
            toy_cohort: Cohort,
            hpo: hpotk.MinimalOntology,
    ):
        # Not Measured and not Observed - 'HP:0006280',  # Chronic pancreatitis
        patient = find_patient('HetSingleVar', toy_cohort)
        term_id = hpotk.TermId.from_curie('HP:0006280')
        predicate = PropagatingPhenotypePredicate(hpo=hpo, query=term_id)
        actual = predicate.test(patient)

        assert actual is None


class TestDiseasePresencePredicate:

    @pytest.mark.parametrize(
        'patient_id, patient_category',
        [
            ('HetSingleVar', PatientCategories.YES),
            ('HomoVar', PatientCategories.NO),
        ])
    def test_disease_predicate(
        self,
        patient_id: str,
        patient_category: PatientCategories,
        toy_cohort: Cohort,

    ):
        patient = find_patient(patient_id, toy_cohort)
        disease_id = hpotk.TermId.from_curie("OMIM:148050")
        predicate = DiseasePresencePredicate(disease_id)
        actual = predicate.test(patient)
        assert actual.phenotype == disease_id
        assert actual.category == patient_category


@pytest.mark.parametrize('patient_id, variant_effect, expected_result',
                         (['HetSingleVar', VariantEffect.FRAMESHIFT_VARIANT, PatientCategories.YES],
                          ['HetSingleVar', VariantEffect.MISSENSE_VARIANT, PatientCategories.NO],
                          ['HetDoubleVar1', VariantEffect.STOP_GAINED, PatientCategories.YES],
                          ['HomoVar', VariantEffect.FRAMESHIFT_VARIANT, PatientCategories.YES],
                          ['LargeCNV', VariantEffect.STOP_LOST, PatientCategories.YES],
                          ['LargeCNV', VariantEffect.FEATURE_TRUNCATION, PatientCategories.YES]))
def test_VariantEffectPredicate(patient_id: str,
                                variant_effect: VariantEffect,
                                expected_result: PatientCategory,
                                toy_cohort: Cohort):
    patient = find_patient(patient_id, toy_cohort)
    predicate = boolean_predicate(VariantPredicates.variant_effect(effect=variant_effect, tx_id='NM_013275.6'))
    result = predicate.test(patient)
    assert result.category == expected_result


@pytest.mark.parametrize('patient_id, variant, hasVarResult',
                         (['HetSingleVar', '16_89279850_89279850_G_GC', PatientCategories.YES],
                          # the `HetSingleVar` patient does NOT have the variant.
                          ['HetSingleVar', '16_89279708_89279725_AGTGTTCGGGGCGGGGCC_A', PatientCategories.NO],
                          # but `HetDoubleVar2` below DOES have the variant.
                          ['HetDoubleVar2', '16_89279708_89279725_AGTGTTCGGGGCGGGGCC_A', PatientCategories.YES],
                          ['HetDoubleVar1', '16_89284601_89284602_GG_A', PatientCategories.YES],
                          ['HetDoubleVar1', '16_89280752_89280752_G_T', PatientCategories.YES],
                          # the `HomoVar` patient does NOT have the variant
                          ['HomoVar', '16_89280752_89280752_G_T', PatientCategories.NO],
                          ['HomoVar', '16_89279458_89279459_TG_T', PatientCategories.YES],
                          ['LargeCNV', '16_89190071_89439815_DEL', PatientCategories.YES]))
def test_VariantKeyPredicate(patient_id, variant, hasVarResult, toy_cohort):
    predicate = boolean_predicate(VariantPredicates.variant_key(key=variant))
    patient = find_patient(patient_id, toy_cohort)
    result = predicate.test(patient)
    assert result.category == hasVarResult


@pytest.mark.parametrize('patient_id, exon, hasVarResult',
                         (['HetSingleVar', 9, PatientCategories.YES],
                          ['HetSingleVar', 13, PatientCategories.NO],
                          ['HetDoubleVar1', 9, PatientCategories.YES],
                          ['HetDoubleVar2', 10, PatientCategories.YES],
                          ['HetDoubleVar2', 9, PatientCategories.YES],
                          ['HomoVar', 10, PatientCategories.NO],
                          ['HomoVar', 9, PatientCategories.YES],
                          ['LargeCNV', 1, PatientCategories.NO],
                          ['LargeCNV', 13, PatientCategories.YES]))
def test_ExonPredicate(patient_id, exon, hasVarResult, toy_cohort):
    patient = find_patient(patient_id, toy_cohort)
    predicate = VariantPredicates.exon(exon=exon, tx_id='NM_013275.6')
    predicate = boolean_predicate(predicate)
    result = predicate.test(patient)
    assert result.category == hasVarResult


@pytest.fixture(scope='module')
def protein_predicates(protein_test_service: ProteinMetadataService) -> ProteinPredicates:
    return ProteinPredicates(protein_metadata_service=protein_test_service)


@pytest.mark.parametrize('patient_id, feature_type, hasVarResult',
                        (['HetDoubleVar2', FeatureType.REGION, PatientCategories.YES],
                        ['HetDoubleVar2', FeatureType.REPEAT, PatientCategories.NO],
                        ['HetSingleVar', FeatureType.REGION, PatientCategories.YES],
                        ['HomoVar', FeatureType.REGION, PatientCategories.YES],
                        ['HetDoubleVar1', FeatureType.REPEAT, PatientCategories.NO]))
## TODO Why do CNV not show as affecting a feature?
##['LargeCNV', FeatureType.REGION , HETEROZYGOUS]))
def test_ProteinFeatureTypePredicate(patient_id, feature_type, hasVarResult, toy_cohort, protein_predicates):
    patient = find_patient(patient_id, toy_cohort)
    predicate = boolean_predicate(protein_predicates.protein_feature_type(feature_type=feature_type, tx_id='NM_013275.6'))
    result = predicate.test(patient)
    assert result.category == hasVarResult


@pytest.mark.parametrize('patient_id, feature, hasVarResult',
                         (['HetDoubleVar2', 'Disordered', PatientCategories.YES],
                          ['HetDoubleVar2', 'BadFeature', PatientCategories.NO],
                          ['HetSingleVar', 'Disordered', PatientCategories.YES],
                          ['HomoVar', 'Disordered', PatientCategories.YES],
                          ['HetDoubleVar1', 'Disordered', PatientCategories.YES]))
def test_ProteinFeaturePredicate(patient_id, feature, hasVarResult, toy_cohort, protein_predicates):
    predicate = boolean_predicate(protein_predicates.protein_feature(feature_id=feature, tx_id='NM_013275.6'))
    patient = find_patient(patient_id, toy_cohort)
    result = predicate.test(patient)
    assert result.category == hasVarResult

@pytest.mark.parametrize('patient_id, region, hasVarResult',
                         (['HetDoubleVar2', Region(2000, 2500), PatientCategories.YES],
                          ['HetDoubleVar2', Region(1000, 1500), PatientCategories.NO],
                          ['HomoVar', Region(2000, 2500), PatientCategories.YES],
                          ['HetSingleVar', Region(2000, 2500), PatientCategories.YES],
                          ['HetDoubleVar1', Region(600, 650), PatientCategories.YES]))
def test_ProteinRegionPredicate(patient_id, region, hasVarResult, toy_cohort):
    predicate = boolean_predicate(VariantPredicates.region(region=region, tx_id='NM_013275.6'))
    patient = find_patient(patient_id, toy_cohort)
    result = predicate.test(patient)
    assert result.category == hasVarResult
