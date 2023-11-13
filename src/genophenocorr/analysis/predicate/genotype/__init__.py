from ._geno import ProtFeaturePredicate, ProtFeatureTypePredicate
from ._geno import VariantEffectPredicate, VariantPredicate, ExonPredicate
from ._geno import HOMOZYGOUS, HETEROZYGOUS, NO_VARIANT

__all__ = [
    'VariantEffectPredicate', 'VariantPredicate', 'ExonPredicate',
    'ProtFeaturePredicate', 'ProtFeatureTypePredicate',
    'HOMOZYGOUS', 'HETEROZYGOUS', 'NO_VARIANT'
]