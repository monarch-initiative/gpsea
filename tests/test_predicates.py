import typing

import hpotk
import pytest

from genophenocorr.analysis.predicate import BooleanPredicate, PatientCategory
from genophenocorr.analysis.predicate.genotype import *
from genophenocorr.analysis.predicate.phenotype import PropagatingPhenotypePredicate
from genophenocorr.model import Cohort, Patient, FeatureType, VariantEffect

from .fixtures import toy_cohort, protein_test_service


def find_patient(pat_id: str, cohort: Cohort) -> typing.Optional[Patient]:
    for pat in cohort.all_patients:
        if pat.patient_id == pat_id:
            return pat


@pytest.mark.parametrize('query, patient_id, expected',
                        # Patient "HetSingleVar" has Phenotypes:
                        # Measured and Observed - 'HP:0001166;HP:0002266',  # Arachnodactyly;Focal clonic seizure
                        # Measured but not Observed - 'HP:0001257',  # Spasticity
                        # Not Measured and not Observed - 'HP:0006280',  # Chronic pancreatitis
                         [
                             # Test exact match
                             ('HP:0001166',  # Arachnodactyly
                             'HetSingleVar',
                              PropagatingPhenotypePredicate.PRESENT),
                             # Test inferred annotations
                             ('HP:0001250',  # Seizure
                              'HetSingleVar',
                              PropagatingPhenotypePredicate.PRESENT),
                             # Test excluded feature
                             ('HP:0001257',  # Spasticity
                              'HetSingleVar',
                              PropagatingPhenotypePredicate.EXCLUDED),
                            ('HP:0006280',  # Chronic pancreatitis
                              'HetSingleVar',
                             PropagatingPhenotypePredicate.NOT_MEASURED),
                         ])
def test_HPOPresentPredicate(toy_cohort: Cohort,
                             toy_hpo: hpotk.Ontology,
                             query: str,
                             patient_id: str,
                             expected: PatientCategory):
    patient = find_patient(patient_id, toy_cohort)
    predicate = PropagatingPhenotypePredicate(hpo=toy_hpo, phenotypic_feature=hpotk.TermId.from_curie(query))
    actual = predicate.test(patient)
    assert actual == expected


@pytest.mark.parametrize('patient_id, variant_effect, expected_result',
                         (['HetSingleVar', VariantEffect.FRAMESHIFT_VARIANT, BooleanPredicate.YES],
                        ['HetSingleVar', VariantEffect.MISSENSE_VARIANT, BooleanPredicate.NO],
                        ['HetDoubleVar1', VariantEffect.STOP_GAINED, BooleanPredicate.YES],
                        ['HomoVar', VariantEffect.FRAMESHIFT_VARIANT, BooleanPredicate.YES],
                        ['LargeCNV', VariantEffect.STOP_LOST, BooleanPredicate.YES],
                        ['LargeCNV', VariantEffect.FEATURE_TRUNCATION, BooleanPredicate.YES]))
def test_VariantEffectPredicate(patient_id: str,
                                variant_effect: VariantEffect,
                                expected_result: PatientCategory,
                                toy_cohort: Cohort):
    patient = find_patient(patient_id, toy_cohort)
    predicate = VariantEffectPredicate('NM_013275.6', effect=variant_effect)
    result = predicate.test(patient)
    assert result == expected_result


@pytest.mark.parametrize('patient_id, variant, hasVarResult',
                         (['HetSingleVar', '16_89279850_89279850_G_GC', BooleanPredicate.YES],
                        # the `HetSingleVar` patient does NOT have the variant.
                        ['HetSingleVar', '16_89279708_89279725_AGTGTTCGGGGCGGGGCC_A', BooleanPredicate.NO],
                        # but `HetDoubleVar2` below DOES have the variant.
                        ['HetDoubleVar2', '16_89279708_89279725_AGTGTTCGGGGCGGGGCC_A', BooleanPredicate.YES],
                        ['HetDoubleVar1', '16_89284601_89284602_GG_A', BooleanPredicate.YES],
                        ['HetDoubleVar1', '16_89280752_89280752_G_T', BooleanPredicate.YES],
                        # the `HomoVar` patient does NOT have the variant
                        ['HomoVar', '16_89280752_89280752_G_T', BooleanPredicate.NO],
                        ['HomoVar', '16_89279458_89279459_TG_T', BooleanPredicate.YES],
                        ['LargeCNV', '16_89190071_89439815_DEL', BooleanPredicate.YES]))
def test_VariantPredicate(patient_id, variant, hasVarResult, toy_cohort):
    predicate = VariantPredicate(variant_key=variant)
    patient = find_patient(patient_id, toy_cohort)
    result = predicate.test(patient)
    assert result == hasVarResult


@pytest.mark.parametrize('patient_id, exon, hasVarResult',
                         (['HetSingleVar', 9, BooleanPredicate.YES],
                        ['HetSingleVar', 13, BooleanPredicate.NO],
                        ['HetDoubleVar1', 9, BooleanPredicate.YES],
                        ['HetDoubleVar2', 10, BooleanPredicate.YES],
                        ['HetDoubleVar2', 9, BooleanPredicate.YES],
                        ['HomoVar', 10, BooleanPredicate.NO],
                        ['HomoVar', 9, BooleanPredicate.YES],
                        ['LargeCNV', 1, BooleanPredicate.NO],
                        ['LargeCNV', 13, BooleanPredicate.YES]))
def test_ExonPredicate(patient_id, exon, hasVarResult, toy_cohort):
    patient = find_patient(patient_id, toy_cohort)
    predicate = ExonPredicate('NM_013275.6', exon_number=exon)
    result = predicate.test(patient)
    assert result == hasVarResult


@pytest.mark.parametrize('patient_id, feature_type, hasVarResult',
                         (['HetDoubleVar2', FeatureType.REGION, BooleanPredicate.YES],
                        ['HetDoubleVar2', FeatureType.REPEAT, BooleanPredicate.NO],
                        ['HetSingleVar', FeatureType.REGION, BooleanPredicate.YES],
                        ['HomoVar', FeatureType.REGION, BooleanPredicate.YES],
                        ['HetDoubleVar1', FeatureType.REPEAT, BooleanPredicate.NO]))
                        ## TODO Why do CNV not show as affecting a feature?
                        ##['LargeCNV', FeatureType.REGION , HETEROZYGOUS]))
def test_ProteinFeatureTypePredicate(patient_id, feature_type, hasVarResult, toy_cohort, protein_test_service):
    patient = find_patient(patient_id, toy_cohort)
    predicate = ProtFeatureTypePredicate('NM_013275.6', feature_type=feature_type, protein_service=protein_test_service)
    result = predicate.test(patient)
    assert result == hasVarResult


@pytest.mark.parametrize('patient_id, feature, hasVarResult',
                         (['HetDoubleVar2', 'Disordered', BooleanPredicate.YES],
                        ['HetDoubleVar2', 'BadFeature', BooleanPredicate.NO],
                        ['HetSingleVar', 'Disordered', BooleanPredicate.YES],
                        ['HomoVar', 'Disordered', BooleanPredicate.YES],
                        ['HetDoubleVar1', 'Disordered', BooleanPredicate.YES]))
def test_ProteinFeaturePredicate(patient_id, feature, hasVarResult, toy_cohort, protein_test_service):
    predicate = ProtFeaturePredicate('NM_013275.6', protein_feature_name=feature, protein_service=protein_test_service)
    patient = find_patient(patient_id, toy_cohort)
    result = predicate.test(patient)
    assert result == hasVarResult


class TestGenericPredicateMethods:

    def test_get_categories(self):
        assert len(PropagatingPhenotypePredicate.get_categories()) == PropagatingPhenotypePredicate.n_categories()
        assert len(BooleanPredicate.get_categories()) == BooleanPredicate.n_categories()
