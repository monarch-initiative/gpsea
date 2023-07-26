from genophenocorr.predicate import *
from genophenocorr.phenotype import PhenotypeCreator
from genophenocorr.constants import VariantEffect
from genophenocorr.protein import FeatureType
import hpotk
import pickle
import os
import pytest

def CreateSmallCohort():
    with open(os.path.join(os.getcwd(), 'tests/samplePatients/testingCohort.pickle'), 'rb') as fh:
        return pickle.load(fh)    

@pytest.fixture
def PatientList():
    patList = CreateSmallCohort().all_patients
    patDict = {}
    for pat in patList:
        patDict[pat.patient_id] = pat
    return patDict

@pytest.fixture
def PhenotypeTest():
    hpo = hpotk.ontology.load.obographs.load_ontology(os.path.join(os.getcwd(), 'tests/testingDefaults/hp.json'))
    validators = [
        hpotk.validate.AnnotationPropagationValidator(hpo),
        hpotk.validate.ObsoleteTermIdsValidator(hpo),
        hpotk.validate.PhenotypicAbnormalityValidator(hpo)
    ]
    phenoCreator = PhenotypeCreator(hpo, hpotk.validate.ValidationRunner(validators))
    return phenoCreator.create_phenotype([('HP:0000729', True)])

@pytest.fixture
def HPOPresentTest():
    return HPOPresentPredicate()

@pytest.mark.parametrize('patient_id, patientCategory', 
                        (['KBG7', PatientCategory.OBSERVED],
                        ['KBG6', PatientCategory.NOTOBSERVED],
                        ['KBG6', PatientCategory.NOTINCLUDED],
                        ['KBG4', PatientCategory.NOTMEASURED],
                        ['KBG4', PatientCategory.NOTINCLUDED]))
def test_HPOPresentPredicate(patient_id, patientCategory, HPOPresentTest, PatientList, PhenotypeTest):
    result = HPOPresentTest.test(PatientList[patient_id], PhenotypeTest[0])
    assert result in patientCategory


@pytest.fixture
def VariantEffectTest():
    return VariantEffectPredicate('NM_013275.6')

@pytest.mark.parametrize('patient_id, variantEffect, hasVarResult', 
                        (['KBG6', VariantEffect.FRAMESHIFT_VARIANT, HasVariantResults.HETEROVARIANT],
                        ['KBG6', VariantEffect.MISSENSE_VARIANT, HasVariantResults.NOVARIANT],
                        ['KBG7', VariantEffect.STOP_GAINED, HasVariantResults.HETEROVARIANT],
                        ['KBG4', VariantEffect.FRAMESHIFT_VARIANT, HasVariantResults.HOMOVARIANT],
                        ['KBG4', VariantEffect.FRAMESHIFT_VARIANT, HasVariantResults.DOMINANTVARIANTS],
                        ['KBG2', VariantEffect.STOP_LOST, HasVariantResults.DOMINANTVARIANTS],
                        ['KBG2', VariantEffect.FEATURE_TRUNCATION, HasVariantResults.DOMINANTVARIANTS]))
def test_VariantEffectPredicate(patient_id, variantEffect, hasVarResult, VariantEffectTest, PatientList):
    result = VariantEffectTest.test(PatientList[patient_id], variantEffect)
    assert result in hasVarResult


@pytest.fixture
def VariantTest():
    return VariantPredicate('NM_013275.6')

@pytest.mark.parametrize('patient_id, variant, hasVarResult', 
                        (['KBG6', '16_89279851_-/C', HasVariantResults.HETEROVARIANT],
                        ['KBG6', '16_89279708_AGTGTTCGGGGCGGGGCC/A', HasVariantResults.NOVARIANT],
                        ['KBG7', '16_89284601_GG/A', HasVariantResults.HETEROVARIANT],
                        ['KBG7', '16_89280752_G/T', HasVariantResults.HETEROVARIANT],
                        ['KBG4', '16_89280752_G/T', HasVariantResults.NOVARIANT],
                        ['KBG4', '16_89279458_TG/T', HasVariantResults.HOMOVARIANT],
                        ['KBG4', '16_89279458_TG/T', HasVariantResults.DOMINANTVARIANTS],
                        ['KBG2', '16_89190071_deletion', HasVariantResults.DOMINANTVARIANTS]))
def test_VariantPredicate(patient_id, variant, hasVarResult, VariantTest, PatientList):
    result = VariantTest.test(PatientList[patient_id], variant)
    assert result in hasVarResult


@pytest.fixture
def ExonTest():
    return ExonPredicate('NM_013275.6')

@pytest.mark.parametrize('patient_id, exon, hasVarResult', 
                        (['KBG6', 9, HasVariantResults.HETEROVARIANT],
                        ['KBG6', 13, HasVariantResults.NOVARIANT],
                        ['KBG7', 9, HasVariantResults.HOMOVARIANT],
                        ['KBG46', 10, HasVariantResults.HETEROVARIANT],
                        ['KBG46', 9, HasVariantResults.HETEROVARIANT],
                        ['KBG4', 10, HasVariantResults.NOVARIANT],
                        ['KBG4', 9, HasVariantResults.DOMINANTVARIANTS],
                        ['KBG4', 9, HasVariantResults.HOMOVARIANT],
                        ['KBG2', 1, HasVariantResults.NOVARIANT],
                        ['KBG2', 13, HasVariantResults.DOMINANTVARIANTS]))
def test_ExonPredicate(patient_id, exon, hasVarResult, ExonTest, PatientList):
    result = ExonTest.test(PatientList[patient_id], exon)
    assert result in hasVarResult


@pytest.fixture
def ProteinFeatureTypeTest():
    return ProtFeatureTypePredicate('NM_013275.6')

@pytest.mark.parametrize('patient_id, featureType, hasVarResult', 
                        (['KBG46', FeatureType.REGION, HasVariantResults.HETEROVARIANT],
                        ['KBG46', FeatureType.REPEAT, HasVariantResults.NOVARIANT],
                        ['KBG6', FeatureType.DOMAIN, HasVariantResults.NOVARIANT],
                        ['KBG4', FeatureType.REGION, HasVariantResults.HOMOVARIANT],
                        ['KBG7', FeatureType.REPEAT, HasVariantResults.NOVARIANT]))
                        ## TODO Why do CNV not show as affecting a feature? 
                        ##['KBG2', FeatureType.REGION , HasVariantResults.DOMINANTVARIANTS]))
def test_ProteinFeatureTypePredicate(patient_id, featureType, hasVarResult, ProteinFeatureTypeTest, PatientList):
    result = ProteinFeatureTypeTest.test(PatientList[patient_id], featureType)
    assert result in hasVarResult


@pytest.fixture
def ProteinFeatureTest():
    return ProtFeaturePredicate('NM_013275.6')

@pytest.mark.parametrize('patient_id, feature, hasVarResult', 
                        (['KBG46', 'Disordered', HasVariantResults.HETEROVARIANT],
                        ['KBG46', 'BadFeature', HasVariantResults.NOVARIANT],
                        ['KBG6', 'Disordered', HasVariantResults.NOVARIANT],
                        ['KBG4', 'Disordered', HasVariantResults.HOMOVARIANT],
                        ['KBG7', 'Disordered', HasVariantResults.NOVARIANT]))
def test_ProteinFeaturePredicate(patient_id, feature, hasVarResult, ProteinFeatureTest, PatientList):
    result = ProteinFeatureTest.test(PatientList[patient_id], feature)
    assert result in hasVarResult

