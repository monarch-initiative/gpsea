import typing

import hpotk
import pytest

from genophenocorr.analysis.predicate import *
from genophenocorr.model import Cohort, Patient, FeatureType, VariantEffect
from .fixtures import toy_hpo, toy_cohort


def find_patient(pat_id, cohort) -> typing.Optional[Patient]:
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
                        (['HetSingleVar', VariantEffect.FRAMESHIFT_VARIANT, BooleanPredicate.TRUE],
                        ['HetSingleVar', VariantEffect.MISSENSE_VARIANT, BooleanPredicate.FALSE],
                        ['HetDoubleVar1', VariantEffect.STOP_GAINED, BooleanPredicate.TRUE],
                        ['HomoVar', VariantEffect.FRAMESHIFT_VARIANT, BooleanPredicate.TRUE],
                        ['LargeCNV', VariantEffect.STOP_LOST, BooleanPredicate.TRUE],
                        ['LargeCNV', VariantEffect.FEATURE_TRUNCATION, BooleanPredicate.TRUE]))
def test_VariantEffectPredicate(patient_id: str,
                                variant_effect: VariantEffect,
                                expected_result: PatientCategory,
                                toy_cohort: Cohort):
    patient = find_patient(patient_id, toy_cohort)
    predicate = VariantEffectPredicate('NM_013275.6', effect=variant_effect)
    result = predicate.test(patient)
    assert result == expected_result


@pytest.mark.parametrize('patient_id, variant, hasVarResult',
                        (['HetSingleVar', '16_89279850_89279850_G_GC', BooleanPredicate.TRUE],
                        # the `HetSingleVar` patient does NOT have the variant.
                        ['HetSingleVar', '16_89279708_89279725_AGTGTTCGGGGCGGGGCC_A', BooleanPredicate.FALSE],
                        # but `HetDoubleVar2` below DOES have the variant.
                        ['HetDoubleVar2', '16_89279708_89279725_AGTGTTCGGGGCGGGGCC_A', BooleanPredicate.TRUE],
                        ['HetDoubleVar1', '16_89284601_89284602_GG_A', BooleanPredicate.TRUE],
                        ['HetDoubleVar1', '16_89280752_89280752_G_T', BooleanPredicate.TRUE],
                        # the `HomoVar` patient does NOT have the variant
                        ['HomoVar', '16_89280752_89280752_G_T', BooleanPredicate.FALSE],
                        ['HomoVar', '16_89279458_89279459_TG_T', BooleanPredicate.TRUE],
                        ['LargeCNV', '16_89190071_89439815_DEL', BooleanPredicate.TRUE]))
def test_VariantPredicate(patient_id, variant, hasVarResult, toy_cohort):
    predicate = VariantPredicate(variant_key=variant)
    patient = find_patient(patient_id, toy_cohort)
    result = predicate.test(patient)
    assert result == hasVarResult


@pytest.mark.parametrize('patient_id, exon, hasVarResult',
                        (['HetSingleVar', 9, BooleanPredicate.TRUE],
                        ['HetSingleVar', 13, BooleanPredicate.FALSE],
                        ['HetDoubleVar1', 9, BooleanPredicate.TRUE],
                        ['HetDoubleVar2', 10, BooleanPredicate.TRUE],
                        ['HetDoubleVar2', 9, BooleanPredicate.TRUE],
                        ['HomoVar', 10, BooleanPredicate.FALSE],
                        ['HomoVar', 9, BooleanPredicate.TRUE],
                        ['LargeCNV', 1, BooleanPredicate.FALSE],
                        ['LargeCNV', 13, BooleanPredicate.TRUE]))
def test_ExonPredicate(patient_id, exon, hasVarResult, toy_cohort):
    patient = find_patient(patient_id, toy_cohort)
    predicate = ExonPredicate('NM_013275.6', exon_number=exon)
    result = predicate.test(patient)
    assert result == hasVarResult


@pytest.mark.parametrize('patient_id, feature_type, hasVarResult',
                        (['HetDoubleVar2', FeatureType.REGION, BooleanPredicate.TRUE],
                        ['HetDoubleVar2', FeatureType.REPEAT, BooleanPredicate.FALSE],
                        ['HetSingleVar', FeatureType.REGION, BooleanPredicate.TRUE],
                        ['HomoVar', FeatureType.REGION, BooleanPredicate.TRUE],
                        ['HetDoubleVar1', FeatureType.REPEAT, BooleanPredicate.FALSE]))
                        ## TODO Why do CNV not show as affecting a feature?
                        ##['LargeCNV', FeatureType.REGION , HETEROZYGOUS]))
def test_ProteinFeatureTypePredicate(patient_id, feature_type, hasVarResult, toy_cohort):
    patient = find_patient(patient_id, toy_cohort)
    predicate = ProtFeatureTypePredicate('NM_013275.6', feature_type=feature_type)
    result = predicate.test(patient)
    assert result == hasVarResult


@pytest.mark.parametrize('patient_id, feature, hasVarResult',
                        (['HetDoubleVar2', 'Disordered', BooleanPredicate.TRUE],
                        ['HetDoubleVar2', 'BadFeature', BooleanPredicate.FALSE],
                        ['HetSingleVar', 'Disordered', BooleanPredicate.TRUE],
                        ['HomoVar', 'Disordered', BooleanPredicate.TRUE],
                        ['HetDoubleVar1', 'Disordered', BooleanPredicate.TRUE]))
def test_ProteinFeaturePredicate(patient_id, feature, hasVarResult, toy_cohort):
    predicate = ProtFeaturePredicate('NM_013275.6', protein_feature_name=feature)
    patient = find_patient(patient_id, toy_cohort)
    result = predicate.test(patient)
    assert result == hasVarResult

