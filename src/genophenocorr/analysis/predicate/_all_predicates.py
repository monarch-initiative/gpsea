import typing

import hpotk

from genophenocorr.model import Patient, FeatureType, VariantEffect

from ._api import PatientCategory, PolyPredicate, BooleanPredicate


class HPOPresentPredicate(PolyPredicate):
    """
    `HPOPresentPredicate` tests if the `patient` is annotated with a `query` HPO term.

    The predicate returns the following results:

    * :attr:`HPOPresentPredicate.PRESENT` if the patient is annotated with the `query` term or its descendant
    * :attr:`HPOPresentPredicate.EXCLUDED` presence of the `query` term or its ancestor was specifically
      excluded in the patient
    * :attr:`HPOPresentPredicate.NOT_MEASURED` if the patient is not annotated with the `query` and presence of `query`
      was *not* excluded
    """

    PRESENT = PatientCategory(cat_id=0,
                              name='Present',
                              description="""
                              The sample *is* annotated with the tested phenotype feature `q`.
                              
                              This is either because the sample is annotated with `q` (exact match),
                              or because one of sample's annotations is a descendant `q` (annotation propagation).
                              For instance, we tested for a Seizure and the sample *had* a Clonic seizure 
                              (a descendant of Seizure).
                              """) #: :meta hide-value:

    EXCLUDED = PatientCategory(cat_id=1,
                               name='Excluded',
                               description="""
                               We are particular about the sample *not* having the tested feature `q`.
                               
                               In other words, `q` was *excluded* in the sample or the sample is annotated with an excluded ancestor of `q`.
                               
                               For instance, we tested for a Clonic seizure and the sample did *not* have any Seizure, which implies 
                               *not* Clonic seizure.  
                               """) #: :meta hide-value:

    NOT_MEASURED = PatientCategory(cat_id=2,
                                   name='Not measured',
                                   description="""
                                   We do not know if the sample has or has not the tested feature.
                                   """) #: :meta hide-value:

    def __init__(self, hpo: hpotk.MinimalOntology,
                 phenotypic_feature: hpotk.TermId) -> None:
        self._hpo = hpotk.util.validate_instance(hpo, hpotk.MinimalOntology, 'hpo')
        self._query = hpotk.util.validate_instance(phenotypic_feature, hpotk.TermId, 'phenotypic_feature')

    def categories(self) -> typing.Sequence[PatientCategory]:
        return HPOPresentPredicate.PRESENT, HPOPresentPredicate.EXCLUDED, HPOPresentPredicate.NOT_MEASURED

    def get_question(self) -> str:
        query_label = self._hpo.get_term(self._query).name
        return f'Is \'{query_label}\' present in the patient?'

    def test(self, patient: Patient) -> typing.Optional[PatientCategory]:
        self._check_patient(patient)
        query = hpotk.util.validate_instance(self._query, hpotk.TermId, 'query')

        if len(patient.phenotypes) == 0:
            return None

        for phenotype in patient.phenotypes:
            if phenotype.is_observed:
                if any(query == anc for anc in self._hpo.graph.get_ancestors(phenotype, include_source=True)):
                    return HPOPresentPredicate.PRESENT
            else:
                if any(query == desc for desc in self._hpo.graph.get_descendants(phenotype, include_source=True)):
                    return self.EXCLUDED

        return self.NOT_MEASURED


HETEROZYGOUS = PatientCategory(cat_id=0,
                               name='Heterozygous',
                               description="""
                               This sample has the tested attribute on one allele. 
                               """)  #: :meta hide-value:

HOMOZYGOUS = PatientCategory(cat_id=1,
                             name='Homozygous',
                             description="""
                             This sample has the tested attribute on both alleles. 
                             """)  #: :meta hide-value:

NO_VARIANT = PatientCategory(cat_id=2,
                             name='No Variant',
                             description="""
                             The sample does not have the tested attribute.
                             """)  #: :meta hide-value:



class VariantEffectPredicate(BooleanPredicate):
    """
    `VariantEffectPredicate` tests if the `patient` has at least one variant that is predicted to have
    the functional `effect` on the transcript of interest.

    :param transcript_id: the accession of the transcript of interest.
    :param effect: the tested variant effect.
    """
    
    def __init__(self, transcript_id: str,
                 effect: VariantEffect) -> None:
        self._tx_id = transcript_id
        self._effect = hpotk.util.validate_instance(effect, VariantEffect, 'effect')

    def get_question(self) -> str:
        return (f'Does the patient have a variant with {self._effect.name} '
                f'on transcript `{self._tx_id}`?')

    def test(self, patient: Patient) -> typing.Optional[PatientCategory]:
        self._check_patient(patient)

        if len(patient.variants) == 0:
            return None

        for variant in patient.variants:
            for ann in variant.tx_annotations:
                if ann.transcript_id == self._tx_id:
                    for var_eff in ann.variant_effects:
                        if var_eff == self._effect:
                            return BooleanPredicate.TRUE

        return BooleanPredicate.FALSE


class VariantPredicate(BooleanPredicate):
    """
    `VariantPredicate` tests if the `patient` has ar least one allele of the variant described by the `variant_key`.

    .. note::

      You can get the variant key by calling :class:`genophenocorr.model.VariantCoordinates.variant_key`.

    :param variant_key: a `str` with the variant key.
    """

    def __init__(self, variant_key: str) -> None:
        self._variant_key = hpotk.util.validate_instance(variant_key, str, 'variant_key')

    def get_question(self) -> str:
        return f'Does the patient have >=1 allele of the variant `{self._variant_key}`?'

    def test(self, patient: Patient) -> typing.Optional[PatientCategory]:
        self._check_patient(patient)

        if len(patient.variants) == 0:
            return None

        for variant in patient.variants:
            if variant.variant_coordinates.variant_key == self._variant_key:
                return BooleanPredicate.TRUE

        return BooleanPredicate.FALSE


class ExonPredicate(BooleanPredicate):
    """
    `ExonPredicate` tests if the `patient` has a variant that affects *n*-th exon of the transcript of interest.

    .. warning::

      We use 1-based numbering to number the exons, not the usual 0-based numbering of the computer science.
      Therefore, the first exon of the transcript has ``exon_number==1``, the second exon is ``2``, and so on ...

    .. warning::

      We do not check if the `exon_number` spans beyond the number of exons of the given `transcript_id`!
      Therefore, ``exon_number==10,000`` will effectively return :attr:`BooleanPredicate.FALSE` for *all* patients!!! 😱
      Well, at least the patients of the *Homo sapiens sapiens* taxon...

    :param transcript_id: the accession of the transcript of interest.
    :param exon_number: a positive `int` of the target exon.
    """

    def __init__(self, transcript_id: str,
                 exon_number: int) -> None:
        self._tx_id = transcript_id
        self._exon_number = hpotk.util.validate_instance(exon_number, int, 'exon_number')
        if self._exon_number <= 0:
            raise ValueError(f'`exon_number` must be a positive `int` but got {self._exon_number}')

    def get_question(self) -> str:
        return (f'Does the patient have a variant that affects exon {self._exon_number} '
                f'on transcript `{self._tx_id}`?')

    def test(self, patient: Patient) -> typing.Optional[PatientCategory]:
        self._check_patient(patient)

        if len(patient.variants) == 0:
            return None

        for variant in patient.variants:
            for ann in variant.tx_annotations:
                if ann.transcript_id == self._tx_id:
                    if ann.overlapping_exons is not None:
                        if self._exon_number in ann.overlapping_exons:
                            return BooleanPredicate.TRUE

        return BooleanPredicate.FALSE


class ProtFeatureTypePredicate(BooleanPredicate):
    """
    `ExonPredicate` tests if the `patient` has a variant that affects a :class:`FeatureType`
    in the transcript of interest.

    :param transcript_id: the accession of the transcript of interest.
    :param feature_type: an instance of the target :class:`FeatureType`.
    """

    def __init__(self, transcript_id: str, feature_type: FeatureType) -> None:
        self._tx_id = transcript_id
        self._feature_type = hpotk.util.validate_instance(feature_type, FeatureType, 'feature_type')

    def get_question(self) -> str:
        return (f'Does the patient have a variant that affects any {self._feature_type.name} protein feature '
                f'on transcript `{self._tx_id}`?')

    def test(self, patient: Patient) -> typing.Optional[PatientCategory]:
        self._check_patient(patient)

        if len(patient.variants) == 0:
            return None

        for variant in patient.variants:
            for ann in variant.tx_annotations:
                if ann.transcript_id == self._tx_id:
                    pe_start, pe_end = ann.protein_effect_location
                    if all(item is not None for item in (pe_start, pe_end)):
                        for prot in ann.protein_affected:
                            for feat in prot.protein_features:
                                if feat.feature_type == self._feature_type:
                                    if len(list(range(max(pe_start, feat.info.start), min(pe_end, feat.info.end) + 1))) > 0:
                                        return BooleanPredicate.TRUE

        return BooleanPredicate.FALSE


class ProtFeaturePredicate(BooleanPredicate):
    """
    `ProtFeaturePredicate` tests if the `patient` has a variant that overlaps with a protein feature.

    The predicate needs the name of the protein feature.
    For instance, `EGF-like 2` for `FBN1 <https://www.uniprot.org/uniprotkb/P35555/entry#family_and_domains>`_

    :param transcript_id: the accession of the transcript of interest.
    :param protein_feature_name: a `str` with the name of the protein feature.
    """

    def __init__(self, transcript_id: str, protein_feature_name: str) -> None:
        self._tx_id = transcript_id
        self._pf_name = hpotk.util.validate_instance(protein_feature_name, str, 'protein_feature_name')

    def get_question(self) -> str:
        return (f'Does the patient have a variant that affects a protein feature with name {self._pf_name} '
                f'on transcript `{self._tx_id}`?')

    def test(self, patient: Patient) -> typing.Optional[PatientCategory]:
        self._check_patient(patient)

        if len(patient.variants) == 0:
            return None

        for variant in patient.variants:
            for ann in variant.tx_annotations:
                if ann.transcript_id == self._tx_id:
                    pe_start, pe_end = ann.protein_effect_location
                    if all(item is not None for item in (pe_start, pe_end)):
                        for prot in ann.protein_affected:
                            for feat in prot.protein_features:
                                if feat.info.name == self._pf_name:
                                    if len(list(range(max(pe_start, feat.info.start), min(pe_end, feat.info.end) + 1))) > 0:
                                        return BooleanPredicate.TRUE

        return BooleanPredicate.FALSE
