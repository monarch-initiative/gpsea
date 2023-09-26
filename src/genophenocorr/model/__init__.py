"""
The `genophenocorr.model` package defines data model classes used in genophenocorr. We start with the top-level elements,
such as :class:`Cohort` and :class:`Patient`, and we follow with data classes for phenotype, genotype, and protein info.
"""

from ._cohort import Cohort, Patient
from ._phenotype import Phenotype
from ._protein import FeatureInfo, FeatureType, ProteinFeature, ProteinMetadata
from ._variant import VariantCoordinates, TranscriptAnnotation, Variant

__all__ = [
    'Cohort', 'Patient',
    'Phenotype',
    'Variant', 'TranscriptAnnotation', 'VariantCoordinates',
    'ProteinMetadata', 'ProteinFeature', 'FeatureInfo', 'FeatureType',
]
