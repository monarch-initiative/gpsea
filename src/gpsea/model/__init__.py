"""
The `gpsea.model` package defines data model classes used in GPSEA.
We start with the top-level elements, such as :class:`~gpsea.model.Cohort`
and :class:`~gpsea.model.Patient`,
and we follow with data classes for phenotype, genotype, transcript, and protein info.
"""

from ._base import SampleLabels, Sex
from ._cohort import Cohort, Patient
from ._gt import Genotype, Genotypes, Genotyped
from ._phenotype import Phenotype, Disease, Measurement
from ._protein import FeatureInfo, FeatureType, ProteinFeature, ProteinMetadata
from ._temporal import Age, AgeKind
from ._tx import TranscriptCoordinates
from ._variant import VariantCoordinates, ImpreciseSvInfo, VariantInfo, VariantInfoAware, VariantClass, Variant
from ._variant import TranscriptAnnotation, TranscriptInfoAware, FunctionalAnnotationAware
from ._variant_effects import VariantEffect

__all__ = [
    'Cohort', 'Patient', 'SampleLabels', 'Sex',
    'Age', 'AgeKind',
    'Phenotype', 'Disease', 'Measurement',
    'Variant', 'VariantClass', 'VariantCoordinates', 'ImpreciseSvInfo', 'VariantInfo', 'VariantInfoAware',
    'Genotype', 'Genotypes', 'Genotyped',
    'TranscriptAnnotation', 'VariantEffect', 'TranscriptInfoAware',
    'FunctionalAnnotationAware', 'TranscriptCoordinates',
    'ProteinMetadata', 'ProteinFeature', 'FeatureInfo', 'FeatureType',
]