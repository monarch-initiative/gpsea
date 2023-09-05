from genophenocorr.predicate import *
from genophenocorr.patient import Patient, BasicPatientCreator
from genophenocorr.phenotype import PhenotypeCreator, Phenotype
from genophenocorr.constants import VariantEffect
from genophenocorr.protein import FeatureType, UniprotProteinMetadataService, ProteinAnnotationCache, ProtCachingFunctionalAnnotator
from genophenocorr.variant import VarCachingFunctionalAnnotator, VariantAnnotationCache, VepFunctionalAnnotator
import hpotk
import os
import pytest

@pytest.fixture
def hpo():
    path = os.path.join(os.getcwd(), 'tests/testingDefaults/hp.toy.json')
    return hpotk.ontology.load.obographs.load_ontology(path)

@pytest.fixture
def PhenotypeCreate(hpo: hpotk.MinimalOntology):
    validators = [
        hpotk.validate.AnnotationPropagationValidator(hpo),
        hpotk.validate.ObsoleteTermIdsValidator(hpo),
        hpotk.validate.PhenotypicAbnormalityValidator(hpo)
    ]
    phenoCreator = PhenotypeCreator(hpo, hpotk.validate.ValidationRunner(validators))
    return phenoCreator


def create_phenotypes(hpo: hpotk.MinimalOntology, observed: str, excluded: str):
    phenotypes = []
    for o in observed.split(';'):
        o = hpotk.TermId.from_curie(o)
        name = hpo.get_term(o).name
        phenotypes.append(Phenotype(o, name=name,  observed=True))

    for e in excluded.split(';'):
        e = hpotk.TermId.from_curie(e)
        name = hpo.get_term(e).name
        phenotypes.append(Phenotype(e, name=name, observed=False))

    return phenotypes


@pytest.mark.parametrize('query, observed, excluded, expected',
                         [
                             # Test exact match
                             ('HP:0001166',  # Arachnodactyly
                              'HP:0001166;HP:0002266',  # Arachnodactyly, Focal clonic seizure
                              'HP:0010677',  # not Enuresis nocturna
                              HPOPresentPredicate.OBSERVED),
                             # Test inferred annotations
                             ('HP:0001250',  # Seizure
                              'HP:0001166;HP:0002266',  # Arachnodactyly, Focal clonic seizure
                              'HP:0010677',  # not Enuresis nocturna
                              HPOPresentPredicate.OBSERVED),

                             # Test excluded feature
                             # TODO - query based on the excluded feature.
                             ('HP:0001250',  # Seizure
                              'HP:0001166;HP:0002266',  # Arachnodactyly, Focal clonic seizure
                              'HP:0010677',  # not Enuresis nocturna
                              HPOPresentPredicate.OBSERVED),
                         ])
def test_HPOPresentPredicate(hpo: hpotk.MinimalOntology,
                             query: str,
                             observed: str,
                             excluded: str,
                             expected: PatientCategory):
    predicate = HPOPresentPredicate(hpo=hpo)
    phenotypes = create_phenotypes(hpo, observed, excluded)
    patient = Patient(patient_id='whatever', phenotypes=phenotypes, variants=[], proteins=[])
    actual = predicate.test(patient, query=hpotk.TermId.from_curie(query))
    assert actual == expected


@pytest.fixture
def PatientDict(PhenotypeCreate):
    pat_id_list = {'HetSingleVar':{'hpos': [['HP:0001250', True]], 'vars': [['16', '89279849', '89279851', 'G', 'GC', 'heterozygous']]},
                'HetDoubleVar1':{'hpos': [['HP:0001250', True]], 'vars': [['16', '89284600', '89284602', 'GG', 'A', 'heterozygous'],['16', '89280751', '89280752', 'G', 'T', 'heterozygous']]},
                'HetDoubleVar2':{'hpos': [['HP:0001250', True]], 'vars': [['16', '89275127', '89275128', 'G', 'A', 'heterozygous'], ['16', '89279707', '89279725', 'AGTGTTCGGGGCGGGGCC', 'A', 'heterozygous']]},
                'HomoVar':{'hpos': [['HP:0001250', True]], 'vars': [['16', '89279457', '89279459', 'TG', 'T', 'homozygous']]},
                'LargeCNV':{'hpos': [['HP:0001250', True]], 'vars': [['16', '89190070', '89439815', 'N', '<DEL>', 'heterozygous']]}
                }
    pat_dict = {}
    pheno_creator = PhenotypeCreate
    # TODO[lnrekerle] - clean the method - remove the cachers/services/...
    path = os.path.join(os.getcwd(), 'tests/sampleAnnotations')
    pm = UniprotProteinMetadataService()
    pac = ProteinAnnotationCache(path)
    pfa = ProtCachingFunctionalAnnotator(pac, pm)
    vac = VariantAnnotationCache(path)
    vep = VepFunctionalAnnotator(pfa)
    vfa = VarCachingFunctionalAnnotator(vac, vep)
    pat_creator = BasicPatientCreator(pheno_creator, vfa)
    for pat, val in pat_id_list.items():
        pat_dict[pat] = pat_creator.create_patient(pat, val['hpos'], val['vars'])
    return pat_dict

@pytest.fixture
def VariantEffectTest():
    return VariantEffectPredicate('NM_013275.6')

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
                                VariantEffectTest: VariantEffectPredicate, 
                                PatientDict
                                ):
    result = VariantEffectTest.test(PatientDict[patient_id], variantEffect)
    assert result == expected_result


@pytest.fixture
def VariantTest():
    return VariantPredicate('NM_013275.6')

@pytest.mark.parametrize('patient_id, variant, hasVarResult', 
                        (['HetSingleVar', '16_89279851_-/C', HETEROZYGOUS],
                        ['HetSingleVar', '16_89279708_AGTGTTCGGGGCGGGGCC/A', NO_VARIANT],
                        ['HetDoubleVar1', '16_89284601_GG/A', HETEROZYGOUS],
                        ['HetDoubleVar1', '16_89280752_G/T', HETEROZYGOUS],
                        ['HomoVar', '16_89280752_G/T', NO_VARIANT],
                        ['HomoVar', '16_89279459_G/-', HOMOZYGOUS],
                        ['LargeCNV', '16_89190071_deletion', HETEROZYGOUS]))
def test_VariantPredicate(patient_id, variant, hasVarResult, VariantTest, PatientDict):
    result = VariantTest.test(PatientDict[patient_id], variant)
    assert result == hasVarResult


@pytest.fixture
def ExonTest():
    return ExonPredicate('NM_013275.6')

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
def test_ExonPredicate(patient_id, exon, hasVarResult, ExonTest, PatientDict):
    result = ExonTest.test(PatientDict[patient_id], exon)
    assert result == hasVarResult


@pytest.fixture
def ProteinFeatureTypeTest():
    return ProtFeatureTypePredicate('NM_013275.6')

@pytest.mark.parametrize('patient_id, featureType, hasVarResult', 
                        (['HetDoubleVar2', FeatureType.REGION, HOMOZYGOUS],
                        ['HetDoubleVar2', FeatureType.REPEAT, NO_VARIANT],
                        ['HetSingleVar', FeatureType.REGION, HETEROZYGOUS],
                        ['HomoVar', FeatureType.REGION, HOMOZYGOUS],
                        ['HetDoubleVar1', FeatureType.REPEAT, NO_VARIANT]))
                        ## TODO Why do CNV not show as affecting a feature? 
                        ##['LargeCNV', FeatureType.REGION , HETEROZYGOUS]))
def test_ProteinFeatureTypePredicate(patient_id, featureType, hasVarResult, ProteinFeatureTypeTest, PatientDict):
    result = ProteinFeatureTypeTest.test(PatientDict[patient_id], featureType)   
    assert result == hasVarResult 


@pytest.fixture
def ProteinFeatureTest():
    return ProtFeaturePredicate('NM_013275.6')

@pytest.mark.parametrize('patient_id, feature, hasVarResult', 
                        (['HetDoubleVar2', 'Disordered', HETEROZYGOUS],
                        ['HetDoubleVar2', 'BadFeature', NO_VARIANT],
                        ['HetSingleVar', 'Disordered', HETEROZYGOUS],
                        ['HomoVar', 'Disordered', HOMOZYGOUS],
                        ['HetDoubleVar1', 'Disordered', HETEROZYGOUS]))
def test_ProteinFeaturePredicate(patient_id, feature, hasVarResult, ProteinFeatureTest, PatientDict):
    result = ProteinFeatureTest.test(PatientDict[patient_id], feature)
    assert result == hasVarResult

