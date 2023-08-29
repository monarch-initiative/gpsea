import hpotk

from genophenocorr.predicate import SimplePredicate, PolyPredicate, PatientCategory
from genophenocorr.patient import Patient
#from genophenocorr.cohort import Cohort
from genophenocorr.constants import VariantEffect
from genophenocorr.phenotype import Phenotype
from genophenocorr.variant import VariantCoordinates
from genophenocorr.protein import FeatureType
import typing
from enum import Flag, auto


class HPOPresentPredicate(PolyPredicate):

    OBSERVED = PatientCategory(cat_id=0,
                               name='Observed',
                               description="""
                               The sample *is* annotated with the tested phenotype feature `q`.
                               
                               This is either because the sample is annotated with `q` (exact match),
                               or because one of sample's annotations is a descendant `q` (annotation propagation).
                               For instance, we tested for a Seizure and the sample did *not* have any Seizure.
                               """)

    NOT_OBSERVED = PatientCategory(cat_id=1,
                                   name='Not observed',
                                   description="""
                                   We are particular about the sample *not* having the tested feature `q`.
                                   
                                   In other words, `q` was *excluded* in the sample or the sample is annotated with an excluded ancestor of `q`.
                                   
                                   For instance, we tested for a Clonic seizure and the sample did *not* have any Seizure, which implies 
                                   *not* Clonic seizure.  
                                   """)

    NOT_MEASURED = PatientCategory(cat_id=2,
                                   name='Not measured',
                                   description="""
                                   We do not know if the sample has or has not the tested feature.
                                   """)

    def __init__(self, 
                 hpo: hpotk.MinimalOntology) -> None:
        self._hpo = hpotk.util.validate_instance(hpo, hpotk.MinimalOntology, 'hpo')

    def categories(self) -> typing.Sequence[PatientCategory]:
        return self.OBSERVED, self.NOT_OBSERVED, self.NOT_MEASURED

    def test(self, patient: Patient, query: hpotk.TermId) -> typing.Optional[PatientCategory]:
        if not isinstance(patient, Patient):
            raise ValueError(f"patient must be type Patient but was type {type(patient)}")
        query = hpotk.util.validate_instance(query, hpotk.TermId, 'query')

        if len(patient.phenotypes) == 0:
            return None

        for pheno in patient.phenotypes:
            if pheno.observed is None:
                continue
            if pheno.observed:
                if any((query == ancestor
                        for ancestor in self._hpo.graph.get_ancestors(pheno, include_source=True))):
                    return self.OBSERVED
            else:
                if any((query == descendant
                        for descendant in self._hpo.graph.get_descendants(pheno, include_source=True))):
                    return self.NOT_OBSERVED

        return self.NOT_MEASURED


HETEROZYGOUS = PatientCategory(cat_id=0,
                            name='Heterozygous',
                            description="""
                            This sample has the tested attribute on one allele. 
                            """)

HOMOZYGOUS = PatientCategory(cat_id=1,
                                name='Homozygous',
                                description="""
                                This sample has the tested attribute on both alleles. 
                                """)

NO_VARIANT = PatientCategory(cat_id=2,
                                name='No Variant',
                                description="""
                                The sample does not have the tested attribute.
                                """)



class VariantEffectPredicate(PolyPredicate):
    
    def __init__(self, transcript:str) -> None:
        self._transcript = transcript

    def categories(self) -> typing.Sequence[PatientCategory]:
        return HETEROZYGOUS, HOMOZYGOUS, NO_VARIANT

    def test(self, patient: Patient, variant_effect:VariantEffect) -> typing.Optional[PatientCategory]:
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
                    return HETEROZYGOUS
                elif v.genotype == "homozygous":
                    return HOMOZYGOUS
                else:
                    return HETEROZYGOUS
        elif len(vars) > 1:
            return HOMOZYGOUS
        else:
            return NO_VARIANT

class VariantPredicate(PolyPredicate):

    def __init__(self, transcript:str) -> None:
        self._transcript = transcript

    def categories(self) -> typing.Sequence[PatientCategory]:
        return HETEROZYGOUS, HOMOZYGOUS, NO_VARIANT

    def test(self, patient: Patient, variant: str) -> typing.Optional[PatientCategory]:
        if not isinstance(patient, Patient):
            raise ValueError(f"patient must be type Patient but was type {type(patient)}")
        if not isinstance(variant, str):
            raise ValueError(f"variant must be type string but was type {type(variant)}")
        vars = set()
        for var in patient.variants:
            print(f"{var.variant_string} == {variant}")
            if var.variant_string == variant:
                vars.add(var)
        if len(vars) == 1:
            for v in vars:
                if v.genotype == "heterozygous":
                    return HETEROZYGOUS
                elif v.genotype == "homozygous":
                    return HOMOZYGOUS
                else:
                    return HETEROZYGOUS
        elif len(vars) > 1:
            return HOMOZYGOUS
        else:
            return NO_VARIANT

class ExonPredicate(PolyPredicate):

    def __init__(self, transcript:str) -> None:
        self._transcript = transcript

    def categories(self) -> typing.Sequence[PatientCategory]:
        return HETEROZYGOUS, HOMOZYGOUS, NO_VARIANT

    def test(self, patient: Patient, exon: int) -> typing.Optional[PatientCategory]:
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
                    return HETEROZYGOUS
                elif v.genotype == "homozygous":
                    return HOMOZYGOUS
                else:
                    return HETEROZYGOUS
        elif len(vars) > 1:
            return HOMOZYGOUS
        else:
            return NO_VARIANT

class ProtFeatureTypePredicate(PolyPredicate):

    def __init__(self, transcript:str) -> None:
         self._transcript = transcript

    def categories(self) -> typing.Sequence[PatientCategory]:
        return HETEROZYGOUS, HOMOZYGOUS, NO_VARIANT

    def test(self, patient: Patient, feature:FeatureType) -> typing.Optional[PatientCategory]:
        if not isinstance(patient, Patient):
            raise ValueError(f"patient must be type Patient but was type {type(patient)}")
        if not isinstance(feature, FeatureType):
            raise ValueError(f"feature must be type FeatureType but was type {type(feature)}")
        vars = set()
        for var in patient.variants:
            for trans in var.tx_annotations:
                if trans.transcript_id == self._transcript:
                    if trans.protein_effect_location is not None and trans.protein_effect_location[0] is not None and trans.protein_effect_location[1] is not None:
                        for prot in trans.protein_affected:
                            for feat in prot.protein_features:
                                if feat.feature_type == feature:
                                    if len(list(range(max(trans.protein_effect_location[0], feat.info.start), min(trans.protein_effect_location[1], feat.info.end) + 1))) > 0:
                                        vars.add(var)
        if len(vars) == 1:
            for v in vars:
                if v.genotype == "heterozygous":
                    return HETEROZYGOUS
                elif v.genotype == "homozygous":
                    return HOMOZYGOUS
                else:
                    return HETEROZYGOUS
        elif len(vars) > 1:
            return HOMOZYGOUS
        else:
            return NO_VARIANT

class ProtFeaturePredicate(PolyPredicate):

    def __init__(self, transcript:str) -> None:
         self._transcript = transcript

    def categories(self) -> typing.Sequence[PatientCategory]:
        return HETEROZYGOUS, HOMOZYGOUS, NO_VARIANT

    def test(self, patient: Patient, feature:str) -> typing.Optional[PatientCategory]:
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
                                    if len(list(range(max(trans.protein_effect_location[0], feat.info.start), min(trans.protein_effect_location[1], feat.info.end) + 1))) > 0:
                                        vars.add(var)
        if len(vars) == 1:
            for v in vars:
                if v.genotype == "heterozygous":
                    return HETEROZYGOUS
                elif v.genotype == "homozygous":
                    return HOMOZYGOUS
                else:
                    return HETEROZYGOUS
        elif len(vars) > 1:
            return HOMOZYGOUS
        else:
            return NO_VARIANT