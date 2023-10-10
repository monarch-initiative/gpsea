import typing

import hpotk

from genophenocorr.model import Patient, FeatureType, Genotype, VariantEffect
from ._api import PolyPredicate, PatientCategory


class HPOPresentPredicate(PolyPredicate[hpotk.TermId]):

    OBSERVED = PatientCategory(cat_id=0,
                               name='Observed',
                               description="""
                               The sample *is* annotated with the tested phenotype feature `q`.
                               
                               This is either because the sample is annotated with `q` (exact match),
                               or because one of sample's annotations is a descendant `q` (annotation propagation).
                               For instance, we tested for a Seizure and the sample *had* a Clonic seizure 
                               (a descendant of Seizure).
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



class VariantEffectPredicate(PolyPredicate[VariantEffect]):
    
    def __init__(self, transcript:str) -> None:
        self._transcript = transcript

    def categories(self) -> typing.Sequence[PatientCategory]:
        return HETEROZYGOUS, HOMOZYGOUS, NO_VARIANT

    def test(self, patient: Patient, query: VariantEffect) -> typing.Optional[PatientCategory]:
        if not isinstance(patient, Patient):
            raise ValueError(f"patient must be type Patient but was type {type(patient)}")
        if not isinstance(query, VariantEffect):
            raise ValueError(f"query must be type VariantEffect but was type {type(query)}")
        vars = set()
        for var in patient.variants:
            for trans in var.tx_annotations:
                if trans.transcript_id == self._transcript:
                    for var_eff in trans.variant_effects:
                        if var_eff == query:
                            vars.add(var)
        if len(vars) == 1:
            for v in vars:
                gt = v.genotype_for_sample(patient.patient_id)
                if gt == Genotype.HETEROZYGOUS:
                    return HETEROZYGOUS
                elif gt == Genotype.HOMOZYGOUS_ALTERNATE:
                    return HOMOZYGOUS
                else:
                    # TODO - is this really what we want to return here?
                    return HETEROZYGOUS
        elif len(vars) > 1:
            return HOMOZYGOUS
        else:
            return NO_VARIANT

class VariantPredicate(PolyPredicate[str]):

    def __init__(self, transcript:str) -> None:
        self._transcript = transcript

    def categories(self) -> typing.Sequence[PatientCategory]:
        return HETEROZYGOUS, HOMOZYGOUS, NO_VARIANT

    def test(self, patient: Patient, query: str) -> typing.Optional[PatientCategory]:
        if not isinstance(patient, Patient):
            raise ValueError(f"patient must be type Patient but was type {type(patient)}")
        if not isinstance(query, str):
            raise ValueError(f"query must be type string but was type {type(query)}")
        vars = set()
        for var in patient.variants:
            if var.variant_coordinates.variant_key == query:
                vars.add(var)
        if len(vars) == 1:
            for v in vars:
                gt = v.genotype_for_sample(patient.patient_id)
                if gt == Genotype.HETEROZYGOUS:
                    return HETEROZYGOUS
                elif gt == Genotype.HOMOZYGOUS_ALTERNATE:
                    return HOMOZYGOUS
                else:
                    # TODO - is this really what we want to return here?
                    return HETEROZYGOUS
        elif len(vars) > 1:
            return HOMOZYGOUS
        else:
            return NO_VARIANT

class ExonPredicate(PolyPredicate[int]):

    def __init__(self, transcript:str) -> None:
        self._transcript = transcript

    def categories(self) -> typing.Sequence[PatientCategory]:
        return HETEROZYGOUS, HOMOZYGOUS, NO_VARIANT

    def test(self, patient: Patient, query: int) -> typing.Optional[PatientCategory]:
        if not isinstance(patient, Patient):
            raise ValueError(f"patient must be type Patient but was type {type(patient)}")
        if not isinstance(query, int):
            raise ValueError(f"query must be type integer but was type {type(query)}")
        vars = set()
        for var in patient.variants:
            for trans in var.tx_annotations:
                if trans.transcript_id == self._transcript:
                    if trans.overlapping_exons is not None:
                        if query in trans.overlapping_exons:
                            vars.add(var)
        if len(vars) == 1:
            for v in vars:
                gt = v.genotype_for_sample(patient.patient_id)
                if gt == Genotype.HETEROZYGOUS:
                    return HETEROZYGOUS
                elif gt == Genotype.HOMOZYGOUS_ALTERNATE:
                    return HOMOZYGOUS
                else:
                    # TODO - is this really what we want to return here?
                    return HETEROZYGOUS
        elif len(vars) > 1:
            return HOMOZYGOUS
        else:
            return NO_VARIANT

class ProtFeatureTypePredicate(PolyPredicate[FeatureType]):

    def __init__(self, transcript:str) -> None:
         self._transcript = transcript

    def categories(self) -> typing.Sequence[PatientCategory]:
        return HETEROZYGOUS, HOMOZYGOUS, NO_VARIANT

    def test(self, patient: Patient, query:FeatureType) -> typing.Optional[PatientCategory]:
        if not isinstance(patient, Patient):
            raise ValueError(f"patient must be type Patient but was type {type(patient)}")
        if not isinstance(query, FeatureType):
            raise ValueError(f"query must be type FeatureType but was type {type(query)}")
        vars = set()
        for var in patient.variants:
            for trans in var.tx_annotations:
                if trans.transcript_id == self._transcript:
                    if trans.protein_effect_location is not None and trans.protein_effect_location[0] is not None and trans.protein_effect_location[1] is not None:
                        for prot in trans.protein_affected:
                            for feat in prot.protein_features:
                                if feat.feature_type == query:
                                    if len(list(range(max(trans.protein_effect_location[0], feat.info.start), min(trans.protein_effect_location[1], feat.info.end) + 1))) > 0:
                                        vars.add(var)
        if len(vars) == 1:
            for v in vars:
                gt = v.genotype_for_sample(patient.patient_id)
                if gt == Genotype.HETEROZYGOUS:
                    return HETEROZYGOUS
                elif gt == Genotype.HOMOZYGOUS_ALTERNATE:
                    return HOMOZYGOUS
                else:
                    # TODO - is this really what we want to return here?
                    return HETEROZYGOUS
        elif len(vars) > 1:
            return HOMOZYGOUS
        else:
            return NO_VARIANT

class ProtFeaturePredicate(PolyPredicate[str]):
    """
    Test if the patient has a variant that overlaps with a protein feature.

    The predicate is generic over a string that represents the name of the protein feature.
    For instance, `EGF-like 2` for `FBN1 <https://www.uniprot.org/uniprotkb/P35555/entry#family_and_domains>`_
    """

    def __init__(self, transcript:str) -> None:
         self._transcript = transcript

    def categories(self) -> typing.Sequence[PatientCategory]:
        return HETEROZYGOUS, HOMOZYGOUS, NO_VARIANT

    def test(self, patient: Patient, query: str) -> typing.Optional[PatientCategory]:
        if not isinstance(patient, Patient):
            raise ValueError(f"patient must be type Patient but was type {type(patient)}")
        if not isinstance(query, str):
            raise ValueError(f"query must be type string but was type {type(query)}")
        vars = set()
        for var in patient.variants:
            for trans in var.tx_annotations:
                if trans.transcript_id == self._transcript:
                    if trans.protein_effect_location is not None and trans.protein_effect_location[0] is not None and trans.protein_effect_location[1] is not None:
                        for prot in trans.protein_affected:
                            for feat in prot.protein_features:
                                if feat.info.name == query:
                                    if len(list(range(max(trans.protein_effect_location[0], feat.info.start), min(trans.protein_effect_location[1], feat.info.end) + 1))) > 0:
                                        vars.add(var)
        if len(vars) == 1:
            for v in vars:
                gt = v.genotype_for_sample(patient.patient_id)
                if gt == Genotype.HETEROZYGOUS:
                    return HETEROZYGOUS
                elif gt == Genotype.HOMOZYGOUS_ALTERNATE:
                    return HOMOZYGOUS
                else:
                    # TODO - is this really what we want to return here?
                    return HETEROZYGOUS
        elif len(vars) > 1:
            return HOMOZYGOUS
        else:
            return NO_VARIANT