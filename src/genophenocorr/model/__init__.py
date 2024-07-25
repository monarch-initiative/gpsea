"""
The `genophenocorr.model` package defines data model classes used in genophenocorr. We start with the top-level elements,
such as :class:`Cohort` and :class:`Patient`, and we follow with data classes for phenotype, genotype, transcript,
and protein info.
"""
from . import genome
from ._base import SampleLabels
from ._cohort import Cohort, Patient
from ._gt import Genotype, Genotypes, Genotyped
from ._phenotype import Phenotype, Disease
from ._protein import FeatureInfo, FeatureType, ProteinFeature, ProteinMetadata
from ._tx import TranscriptCoordinates
from ._variant import VariantCoordinates, TranscriptAnnotation, TranscriptInfoAware, Variant
from ._variant_effects import VariantEffect

__all__ = [
    'Cohort', 'Patient', 'SampleLabels',
    'Phenotype', 'Disease',
    'Variant', 'VariantCoordinates', 'Genotype', 'Genotypes', 'Genotyped', 
    'TranscriptAnnotation', 'VariantEffect', 'TranscriptInfoAware', 'TranscriptCoordinates',
    'ProteinMetadata', 'ProteinFeature', 'FeatureInfo', 'FeatureType',
]
