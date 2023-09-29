"""
The `genophenocorr.model` package defines data model classes used in genophenocorr. We start with the top-level elements,
such as :class:`Cohort` and :class:`Patient`, and we follow with data classes for phenotype, genotype, transcript,
and protein info.
"""
from . import genome

from ._cohort import Cohort, Patient
from ._phenotype import Phenotype
from ._protein import FeatureInfo, FeatureType, ProteinFeature, ProteinMetadata
from ._tx import TranscriptCoordinates
from ._variant import VariantCoordinates, TranscriptAnnotation, TranscriptInfoAware, Variant

__all__ = [
    'Cohort', 'Patient',
    'Phenotype',
    'Variant', 'VariantCoordinates', 'TranscriptAnnotation', 'TranscriptInfoAware', 'TranscriptCoordinates',
    'ProteinMetadata', 'ProteinFeature', 'FeatureInfo', 'FeatureType',
]
