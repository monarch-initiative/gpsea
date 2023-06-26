from genophenocorr.predicate import SimplePredicate, PolyPredicate
from genophenocorr.patient import Patient
#from genophenocorr.cohort import Cohort
from genophenocorr.constants import VariantEffect
from genophenocorr.phenotype import Phenotype
from genophenocorr.variant import VariantCoordinates
from genophenocorr.protein import FeatureType
import typing
from collections import namedtuple
import os

PatientCategory = namedtuple('PatientCategory', field_names=['cat_id', 'name'])


class HPOPresentPredicate(PolyPredicate):
    def __init__(self) -> None:
        pass

    def categories(self) -> typing.Sequence[PatientCategory]:
        pass

    def test(self, patient: Patient, hpo:Phenotype) -> typing.Optional[PatientCategory]:
        if not isinstance(hpo, Phenotype):
            raise ValueError(f"hpo must be type Phenotype but was type {type(hpo)}")
        if not isinstance(patient, Patient):
            raise ValueError(f"patient must be type Patient but was type {type(patient)}")
        for pheno in patient.phenotypes:
            if hpo == pheno and pheno.observed == True:
                return PatientCategory('HasHPO', 'Observed')
            elif hpo == pheno and pheno.observed == False:
                return PatientCategory('HasHPO', 'NotObserved')
        return PatientCategory('HasHPO', 'NotMeasured')

class VariantEffectPredicate(SimplePredicate):
    
    def __init__(self, transcript:str) -> None:
        self._transcript = transcript

    def test(self, patient: Patient, variant_effect: VariantEffect) -> bool:
        if not isinstance(patient, Patient):
            raise ValueError(f"patient must be type Patient but was type {type(patient)}")
        if not isinstance(variant_effect, VariantEffect):
            raise ValueError(f"variant_effect must be type VariantEffect but was type {type(variant_effect)}")
        for var in patient.variants:
            for trans in var.tx_annotations:
                if trans.transcript_id == self._transcript:
                    for var_eff in trans.variant_effects:
                        if var_eff == variant_effect.effect_id:
                            return True
        return False

class VariantPredicate(SimplePredicate):

    def __init__(self, transcript:str) -> None:
        self._transcript = transcript

    def test(self, patient: Patient, variant: str) -> bool:
        if not isinstance(patient, Patient):
            raise ValueError(f"patient must be type Patient but was type {type(patient)}")
        if not isinstance(variant, str):
            raise ValueError(f"variant must be type string but was type {type(variant)}")
        for var in patient.variants:
            if var.variant_string == variant:
                return True
        return False

class ExonPredicate(SimplePredicate):

    def __init__(self, transcript:str) -> None:
        self._transcript = transcript

    def test(self, patient: Patient, exon: int) -> bool:
        if not isinstance(patient, Patient):
            raise ValueError(f"patient must be type Patient but was type {type(patient)}")
        if not isinstance(exon, int):
            raise ValueError(f"exon must be type integer but was type {type(exon)}")
        for var in patient.variants:
            for trans in var.tx_annotations:
                if trans.transcript_id == self._transcript:
                    if trans.overlapping_exons is not None:
                        if exon in trans.overlapping_exons:
                            return True
        return False

class ProtFeatureTypePredicate(SimplePredicate):

    def __init__(self, transcript:str) -> None:
         self._transcript = transcript

    def test(self, patient: Patient, feature:FeatureType) -> bool:
        if not isinstance(patient, Patient):
            raise ValueError(f"patient must be type Patient but was type {type(patient)}")
        if not isinstance(feature, FeatureType):
            raise ValueError(f"domain must be type FeatureType but was type {type(feature)}")
        for var in patient.variants:
            for trans in var.tx_annotations:
                if trans.transcript_id == self._transcript:
                    if trans.protein_effect_location is not None and trans.protein_effect_location[0] is not None and trans.protein_effect_location[1] is not None:
                        for prot in trans.protein_affected:
                            for feat in prot.protein_features:
                                if feat.feature_type == feature:
                                    if any(check in range(trans.protein_effect_location[0], trans.protein_effect_location[1]) for check in range(feat.info.start, feat.info.end)):
                                        return True
        return False

class ProtFeaturePredicate(SimplePredicate):

    def __init__(self, transcript:str) -> None:
         self._transcript = transcript

    def test(self, patient: Patient, feature:str) -> bool:
        if not isinstance(patient, Patient):
            raise ValueError(f"patient must be type Patient but was type {type(patient)}")
        if not isinstance(feature, str):
            raise ValueError(f"domain must be type string but was type {type(feature)}")
        for var in patient.variants:
            for trans in var.tx_annotations:
                if trans.transcript_id == self._transcript:
                    if trans.protein_effect_location is not None and trans.protein_effect_location[0] is not None and trans.protein_effect_location[1] is not None:
                        for prot in trans.protein_affected:
                            for feat in prot.protein_features:
                                if feat.info.name == feature:
                                    if any(check in range(trans.protein_effect_location[0], trans.protein_effect_location[1]) for check in range(feat.info.start, feat.info.end)):
                                        return True
        return False