from abc import ABCMeta
from typing import Any
from genophenocorr.model import Variant, VariantEffect, FeatureType
from genophenocorr.model.genome import Region
from ._api import VariantPredicate
from genophenocorr.preprocessing import ProteinMetadataService


class VariantEffectPredicate(VariantPredicate):
    
    def __init__(self, effect: VariantEffect, tx_id: str) -> None:
        self._effect = effect
        self._tx_id = tx_id
        
    def test(self, variant: Variant) -> bool:
        
        tx_anno = variant.get_tx_anno_by_id(self._tx_id)
        for effect in tx_anno.variant_effects:
            if effect == self._effect:
                return True
        return False
    
    
class VariantKeyPredicate(VariantPredicate):
    
    def __init__(self, key: str) -> None:
        self._key = key
        
    def test(self, variant: Variant) -> bool:
        
        if variant.variant_coordinates.variant_key == self._key:
            return True
        return False
    
class VariantGenePredicate(VariantPredicate):
    
    def __init__(self, gene_symbol:str) -> None:
        self._symbol = gene_symbol
        
    def test(self, variant: Variant) -> bool:
        for tx in variant.tx_annotations:
            if tx.gene_id == self._symbol:
                return True
        return False
    
class VariantTranscriptPredicate(VariantPredicate):
    
    def __init__(self, tx_id: str) -> None:
        self._tx_id = tx_id
        
    def test(self, variant: Variant) -> bool:
        for tx in variant.tx_annotations:
            if tx.transcript_id == self._tx_id:
                return True
        return False
    
    
class VariantExonPredicate(VariantPredicate):
    
    def __init__(self, exon: int, tx_id: str) -> None:
        self._exon = exon
        self._tx_id = tx_id
        
    def test(self, variant: Variant) -> bool:
        tx_anno = variant.get_tx_anno_by_id(self._tx_id)
        if self._exon in tx_anno.overlapping_exons:
            return True
        return False
    
class ProteinRegionPredicate(VariantPredicate):
    
    def __init__(self, region: Region, tx_id: str) -> None:
        self._region = region
        self._tx_id = tx_id
        
    def test(self, variant: Variant) -> bool:
        tx_anno = variant.get_tx_anno_by_id(self._tx_id)
        if tx_anno.protein_effect_location.overlaps_with(self._region):
            return True
        return False
    
class ProteinFeatureTypePredicate(VariantPredicate):
    
    def __init__(self, feature_type: FeatureType, tx_id: str, protein_metadata_service: ProteinMetadataService) -> None:
        self._feature_type = feature_type
        self._tx_id = tx_id
        self._prot_service = protein_metadata_service
        
    def test(self, variant: Variant) -> bool:
        tx_anno = variant.get_tx_anno_by_id(self._tx_id)
        protein = self._prot_service.annotate(tx_anno.protein_id)
        for feat in protein.protein_features:
            if feat.feature_type == self._feature_type and tx_anno.protein_effect_location.overlaps_with(feat.info.region):
                return True
        return False
    
class ProteinFeaturePredicate(VariantPredicate):
    
    def __init__(self, feature_name: str, tx_id: str, protein_metadata_service: ProteinMetadataService) -> None:
        self._feature = feature_name
        self._tx_id = tx_id
        self._prot_service = protein_metadata_service
        
    def test(self, variant: Variant) -> bool:
        tx_anno = variant.get_tx_anno_by_id(self._tx_id)
        protein = self._prot_service.annotate(tx_anno.protein_id)
        for feat in protein.protein_features:
            if feat.info.name == self._feature and tx_anno.protein_effect_location.overlaps_with(feat.info.region):
                return True
        return False
    
    