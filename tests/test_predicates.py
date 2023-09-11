import typing

import hpotk
import pytest

from genophenocorr.analysis.predicate import *
from genophenocorr.constants import VariantEffect
from genophenocorr.model import Cohort, Patient, FeatureType
from .test_fixtures import toy_hpo, test_cohort


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
                              HPOPresentPredicate.OBSERVED),
                             # Test inferred annotations
                             ('HP:0001250',  # Seizure
                              'HetSingleVar',
                              HPOPresentPredicate.OBSERVED),
                             # Test excluded feature
                             ('HP:0001257',  # Spasticity
                              'HetSingleVar',
                              HPOPresentPredicate.NOT_OBSERVED),
                            ('HP:0006280',  # Chronic pancreatitis
                              'HetSingleVar',
                              HPOPresentPredicate.NOT_MEASURED),
                         ])
def test_HPOPresentPredicate(test_cohort: Cohort,
                             toy_hpo: hpotk.Ontology,
                             query: str,
                             patient_id: str,
                             expected: PatientCategory):
    patient = find_patient(patient_id, test_cohort)
    predicate = HPOPresentPredicate(hpo=toy_hpo)
    actual = predicate.test(patient, query=hpotk.TermId.from_curie(query))
    assert actual == expected


@pytest.mark.parametrize('patient_id, variantEffect, expected_result',
                        (['HetSingleVar', VariantEffect.FRAMESHIFT_VARIANT, HETEROZYGOUS],
                        ['HetSingleVar', VariantEffect.MISSENSE_VARIANT, NO_VARIANT],
                        ['HetDoubleVar1', VariantEffect.STOP_GAINED, HETEROZYGOUS],
                        ['HomoVar', VariantEffect.FRAMESHIFT_VARIANT, HOMOZYGOUS],
                        ['LargeCNV', VariantEffect.STOP_LOST, HETEROZYGOUS],
                        ['LargeCNV', VariantEffect.FEATURE_TRUNCATION, HETEROZYGOUS]))
def test_VariantEffectPredicate(patient_id: str, 
                                variantEffect: VariantEffect,
                                expected_result: PatientCategory,
                                test_cohort: Cohort):
    patient = find_patient(patient_id, test_cohort)
    predicate = VariantEffectPredicate('NM_013275.6')
    result = predicate.test(patient, variantEffect)
    assert result == expected_result


@pytest.mark.parametrize('patient_id, variant, hasVarResult',
                        (['HetSingleVar', '16_89279851_-/C', HETEROZYGOUS],
                        ['HetSingleVar', '16_89279708_AGTGTTCGGGGCGGGGCC/A', NO_VARIANT],
                        ['HetDoubleVar1', '16_89284601_GG/A', HETEROZYGOUS],
                        ['HetDoubleVar1', '16_89280752_G/T', HETEROZYGOUS],
                        ['HomoVar', '16_89280752_G/T', NO_VARIANT],
                        ['HomoVar', '16_89279458_TG/T', HOMOZYGOUS],
                        ['LargeCNV', '16_89190071_deletion', HETEROZYGOUS]))
def test_VariantPredicate(patient_id, variant, hasVarResult, test_cohort):
    predicate = VariantPredicate('NM_013275.6')
    patient = find_patient(patient_id, test_cohort)
    result = predicate.test(patient, variant)
    assert result == hasVarResult


@pytest.mark.parametrize('patient_id, exon, hasVarResult',
                        (['HetSingleVar', 9, HETEROZYGOUS],
                        ['HetSingleVar', 13, NO_VARIANT],
                        ['HetDoubleVar1', 9, HOMOZYGOUS],
                        ['HetDoubleVar2', 10, HETEROZYGOUS],
                        ['HetDoubleVar2', 9, HETEROZYGOUS],
                        ['HomoVar', 10, NO_VARIANT],
                        ['HomoVar', 9, HOMOZYGOUS],
                        ['LargeCNV', 1, NO_VARIANT],
                        ['LargeCNV', 13, HETEROZYGOUS]))
def test_ExonPredicate(patient_id, exon, hasVarResult, test_cohort):
    patient = find_patient(patient_id, test_cohort)
    predicate = ExonPredicate('NM_013275.6')
    result = predicate.test(patient, exon)
    assert result == hasVarResult


@pytest.mark.parametrize('patient_id, featureType, hasVarResult',
                        (['HetDoubleVar2', FeatureType.REGION, HOMOZYGOUS],
                        ['HetDoubleVar2', FeatureType.REPEAT, NO_VARIANT],
                        ['HetSingleVar', FeatureType.REGION, HETEROZYGOUS],
                        ['HomoVar', FeatureType.REGION, HOMOZYGOUS],
                        ['HetDoubleVar1', FeatureType.REPEAT, NO_VARIANT]))
                        ## TODO Why do CNV not show as affecting a feature?
                        ##['LargeCNV', FeatureType.REGION , HETEROZYGOUS]))
def test_ProteinFeatureTypePredicate(patient_id, featureType, hasVarResult, test_cohort):
    patient = find_patient(patient_id, test_cohort)
    predicate = ProtFeatureTypePredicate('NM_013275.6')
    result = predicate.test(patient, featureType)
    assert result == hasVarResult


@pytest.mark.parametrize('patient_id, feature, hasVarResult',
                        (['HetDoubleVar2', 'Disordered', HETEROZYGOUS],
                        ['HetDoubleVar2', 'BadFeature', NO_VARIANT],
                        ['HetSingleVar', 'Disordered', HETEROZYGOUS],
                        ['HomoVar', 'Disordered', HOMOZYGOUS],
                        ['HetDoubleVar1', 'Disordered', HETEROZYGOUS]))
def test_ProteinFeaturePredicate(patient_id, feature, hasVarResult, test_cohort):
    predicate = ProtFeaturePredicate('NM_013275.6')
    patient = find_patient(patient_id, test_cohort)
    result = predicate.test(patient, feature)
    assert result == hasVarResult

