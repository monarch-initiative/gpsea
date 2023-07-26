from genophenocorr.predicate import SimplePredicate, PolyPredicate
from genophenocorr.patient import Patient
#from genophenocorr.cohort import Cohort
from genophenocorr.constants import VariantEffect
from genophenocorr.phenotype import Phenotype
from genophenocorr.variant import VariantCoordinates
from genophenocorr.protein import FeatureType
import typing
from collections import namedtuple
from enum import Flag, auto
import os

class PatientCategory(Flag):
    OBSERVED = auto()
    NOTOBSERVED = auto()
    NOTMEASURED = auto()
    NOTINCLUDED = NOTOBSERVED | NOTMEASURED


class HPOPresentPredicate(PolyPredicate):
    def __init__(self) -> None:
        self._categories = self.categories()

    def categories(self) -> typing.Sequence[PatientCategory]:
        pass

    def test(self, patient: Patient, hpo:Phenotype) -> typing.Optional[PatientCategory]:
        if not isinstance(hpo, Phenotype):
            raise ValueError(f"hpo must be type Phenotype but was type {type(hpo)}")
        if not isinstance(patient, Patient):
            raise ValueError(f"patient must be type Patient but was type {type(patient)}")
        for pheno in patient.phenotypes:
            if hpo == pheno and pheno.observed == True:
                return PatientCategory.OBSERVED
            elif hpo == pheno and pheno.observed == False:
                return PatientCategory.NOTOBSERVED
        return PatientCategory.NOTMEASURED


class HasVariantResults(Flag):
    NOVARIANT = auto()
    HETEROVARIANT = auto()
    HOMOVARIANT = auto()
    DOMINANTVARIANTS = HETEROVARIANT | HOMOVARIANT



class VariantEffectPredicate(SimplePredicate):
    
    def __init__(self, transcript:str) -> None:
        self._transcript = transcript

    def test(self, patient: Patient, variant_effect:VariantEffect) -> HasVariantResults:
        if not isinstance(patient, Patient):
            raise ValueError(f"patient must be type Patient but was type {type(patient)}")
        if not isinstance(variant_effect, VariantEffect):
            raise ValueError(f"variant_effect must be type VariantEffect but was type {type(variant_effect)}")
        vars = set()
        for var in patient.variants:
            for trans in var.tx_annotations:
                if trans.transcript_id == self._transcript:
                    for var_eff in trans.variant_effects:
                        if var_eff == str(variant_effect):
                            vars.add(var)
        if len(vars) == 1:
            for v in vars:
                if v.genotype == "heterozygous":
                    return HasVariantResults.HETEROVARIANT
                elif v.genotype == "homozygous":
                    return HasVariantResults.HOMOVARIANT
                else:
                    return HasVariantResults.DOMINANTVARIANTS
        elif len(vars) > 1:
            return HasVariantResults.HOMOVARIANT
        else:
            return HasVariantResults.NOVARIANT

class VariantPredicate(SimplePredicate):

    def __init__(self, transcript:str) -> None:
        self._transcript = transcript

    def test(self, patient: Patient, variant: str) -> HasVariantResults:
        if not isinstance(patient, Patient):
            raise ValueError(f"patient must be type Patient but was type {type(patient)}")
        if not isinstance(variant, str):
            raise ValueError(f"variant must be type string but was type {type(variant)}")
        vars = set()
        for var in patient.variants:
            if var.variant_string == variant:
                vars.add(var)
        if len(vars) == 1:
            for v in vars:
                if v.genotype == "heterozygous":
                    return HasVariantResults.HETEROVARIANT
                elif v.genotype == "homozygous":
                    return HasVariantResults.HOMOVARIANT
                else:
                    return HasVariantResults.DOMINANTVARIANTS
        elif len(vars) > 1:
            return HasVariantResults.HOMOVARIANT
        else:
            return HasVariantResults.NOVARIANT

class ExonPredicate(SimplePredicate):

    def __init__(self, transcript:str) -> None:
        self._transcript = transcript

    def test(self, patient: Patient, exon: int) -> HasVariantResults:
        if not isinstance(patient, Patient):
            raise ValueError(f"patient must be type Patient but was type {type(patient)}")
        if not isinstance(exon, int):
            raise ValueError(f"exon must be type integer but was type {type(exon)}")
        vars = set()
        for var in patient.variants:
            for trans in var.tx_annotations:
                if trans.transcript_id == self._transcript:
                    if trans.overlapping_exons is not None:
                        if exon in trans.overlapping_exons:
                            vars.add(var)
        if len(vars) == 1:
            for v in vars:
                if v.genotype == "heterozygous":
                    return HasVariantResults.HETEROVARIANT
                elif v.genotype == "homozygous":
                    return HasVariantResults.HOMOVARIANT
                else:
                    return HasVariantResults.DOMINANTVARIANTS
        elif len(vars) > 1:
            return HasVariantResults.HOMOVARIANT
        else:
            return HasVariantResults.NOVARIANT

class ProtFeatureTypePredicate(SimplePredicate):

    def __init__(self, transcript:str) -> None:
         self._transcript = transcript

    def test(self, patient: Patient, feature:FeatureType) -> HasVariantResults:
        if not isinstance(patient, Patient):
            raise ValueError(f"patient must be type Patient but was type {type(patient)}")
        if not isinstance(feature, FeatureType):
            raise ValueError(f"domain must be type FeatureType but was type {type(feature)}")
        vars = set()
        for var in patient.variants:
            for trans in var.tx_annotations:
                if trans.transcript_id == self._transcript:
                    if trans.protein_effect_location is not None and trans.protein_effect_location[0] is not None and trans.protein_effect_location[1] is not None:
                        for prot in trans.protein_affected:
                            for feat in prot.protein_features:
                                if feat.feature_type == feature:
                                    if any(check in range(trans.protein_effect_location[0], trans.protein_effect_location[1]) for check in range(feat.info.start, feat.info.end)):
                                        vars.add(var)
        if len(vars) == 1:
            for v in vars:
                if v.genotype == "heterozygous":
                    return HasVariantResults.HETEROVARIANT
                elif v.genotype == "homozygous":
                    return HasVariantResults.HOMOVARIANT
                else:
                    return HasVariantResults.DOMINANTVARIANTS
        elif len(vars) > 1:
            return HasVariantResults.HOMOVARIANT
        else:
            return HasVariantResults.NOVARIANT

class ProtFeaturePredicate(SimplePredicate):

    def __init__(self, transcript:str) -> None:
         self._transcript = transcript

    def test(self, patient: Patient, feature:str) -> HasVariantResults:
        if not isinstance(patient, Patient):
            raise ValueError(f"patient must be type Patient but was type {type(patient)}")
        if not isinstance(feature, str):
            raise ValueError(f"domain must be type string but was type {type(feature)}")
        vars = set()
        for var in patient.variants:
            for trans in var.tx_annotations:
                if trans.transcript_id == self._transcript:
                    if trans.protein_effect_location is not None and trans.protein_effect_location[0] is not None and trans.protein_effect_location[1] is not None:
                        for prot in trans.protein_affected:
                            for feat in prot.protein_features:
                                if feat.info.name == feature:
                                    if any(check in range(trans.protein_effect_location[0], trans.protein_effect_location[1]) for check in range(feat.info.start, feat.info.end)):
                                        vars.add(var)
        if len(vars) == 1:
            for v in vars:
                if v.genotype == "heterozygous":
                    return HasVariantResults.HETEROVARIANT
                elif v.genotype == "homozygous":
                    return HasVariantResults.HOMOVARIANT
                else:
                    return HasVariantResults.DOMINANTVARIANTS
        elif len(vars) > 1:
            return HasVariantResults.HOMOVARIANT
        else:
            return HasVariantResults.NOVARIANT