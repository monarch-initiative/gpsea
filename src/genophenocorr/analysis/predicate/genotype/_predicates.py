from abc import ABCMeta
from typing import Any
from genophenocorr.model import Variant, VariantEffect, FeatureType
from genophenocorr.model.genome import Region
from ._api import VariantPredicate
from genophenocorr.preprocessing import ProteinMetadataService


class VariantEffectPredicate(VariantPredicate):
    """
    `VariantEffectPredicate` is a `VariantPredicate` that sets up testing
    for a specific `VariantEffect` on a given transcript ID. 
    
    Args:
        effect (VariantEffect): 
        tx_id (str): 
    """
    # TODO: add __repr__, __str__, __hash__, __eq__
    
    def __init__(self, effect: VariantEffect, tx_id: str) -> None:
        self._effect = effect
        self._tx_id = tx_id
        
    def get_question(self) -> str:
        return f'{self._effect.name} on {self._tx_id}'

    def test(self, variant: Variant) -> bool:
        """Tests if the `Variant` causes the specified `VariantEffect`
        on the given transcript. 

        Args:
            variant (Variant): a 

        Returns:
            bool: _description_
        """
        
        tx_anno = variant.get_tx_anno_by_id(self._tx_id)
        if tx_anno is None:
            return False
        for effect in tx_anno.variant_effects:
            if effect == self._effect:
                return True
        return False
    

class VariantKeyPredicate(VariantPredicate):
    # TODO: add __repr__, __str__, __hash__, __eq__
    
    def __init__(self, key: str) -> None:
        self._key = key
        
    def get_question(self) -> str:
        return f'variant has ID of {self._key}'

    def test(self, variant: Variant) -> bool:
        return self._key == variant.variant_coordinates.variant_key
    
class VariantGenePredicate(VariantPredicate):
    # TODO: add __repr__, __str__, __hash__, __eq__
    
    def __init__(self, gene_symbol:str) -> None:
        self._symbol = gene_symbol

    def get_question(self) -> str:
        return f'variant affects gene {self._symbol}'

    def test(self, variant: Variant) -> bool:
        for tx in variant.tx_annotations:
            if tx.gene_id == self._symbol:
                return True
        return False
    
class VariantTranscriptPredicate(VariantPredicate):
    # TODO: add __repr__, __str__, __hash__, __eq__
    
    def __init__(self, tx_id: str) -> None:
        self._tx_id = tx_id

    def get_question(self) -> str:
        return f'variant affects transcript {self._tx_id}'
        
    def test(self, variant: Variant) -> bool:
        for tx in variant.tx_annotations:
            if tx.transcript_id == self._tx_id:
                return True
        return False
    
    
class VariantExonPredicate(VariantPredicate):
    # TODO: add __repr__, __str__, __hash__, __eq__
    
    def __init__(self, exon: int, tx_id: str) -> None:
        self._exon = exon
        self._tx_id = tx_id
        
    def get_question(self) -> str:
        return f'variant affects exon {self._exon} on {self._tx_id}'

    def test(self, variant: Variant) -> bool:
        tx_anno = variant.get_tx_anno_by_id(self._tx_id)
        if tx_anno is None:
            return False
        
        if tx_anno.overlapping_exons is None:
            return False

        return any(self._exon == exon for exon in tx_anno.overlapping_exons)
    
class ProteinRegionPredicate(VariantPredicate):
    # TODO: add __repr__, __str__, __hash__, __eq__
    
    def __init__(self, region: Region, tx_id: str) -> None:
        self._region = region
        self._tx_id = tx_id
        
    def get_question(self) -> str:
        return f'variant affects aminoacid(s) between {self._region.start} and {self._region.end} on protein encoded by transcript {self._tx_id}'

    def test(self, variant: Variant) -> bool:
        tx_anno = variant.get_tx_anno_by_id(self._tx_id)
        if tx_anno is None:
            return False
        location = tx_anno.protein_effect_location
        if location is None:
            return False
        return location.overlaps_with(self._region)
    
class ProteinFeatureTypePredicate(VariantPredicate):
    # TODO: add __repr__, __str__, __hash__, __eq__
    
    def __init__(self, feature_type: FeatureType, tx_id: str, protein_metadata_service: ProteinMetadataService) -> None:
        self._feature_type = feature_type
        self._tx_id = tx_id
        self._prot_service = protein_metadata_service
        
    def get_question(self) -> str:
        return f'variant affects {self._feature_type.name} feature type on the protein encoded by transcript {self._tx_id}'

    def test(self, variant: Variant) -> bool:
        tx_anno = variant.get_tx_anno_by_id(self._tx_id)
        if tx_anno is None:
            return False
        location = tx_anno.protein_effect_location
        if location is None:
            return False
        
        protein = self._prot_service.annotate(tx_anno.protein_id)
        for feat in protein.protein_features:    
            if feat.feature_type == self._feature_type and location.overlaps_with(feat.info.region):
                return True
        return False
    
class ProteinFeaturePredicate(VariantPredicate):
    # TODO: add __repr__, __str__, __hash__, __eq__
    
    def __init__(self, feature_name: str, tx_id: str, protein_metadata_service: ProteinMetadataService) -> None:
        self._feature = feature_name
        self._tx_id = tx_id
        self._prot_service = protein_metadata_service

    def get_question(self) -> str:
        return f'Variant that affects {self._feature} feature on the protein encoded by transcript {self._tx_id}'
        
    def test(self, variant: Variant) -> bool:
        tx_anno = variant.get_tx_anno_by_id(self._tx_id)
        if tx_anno is None:
            return False
        
        location = tx_anno.protein_effect_location
        if location is None:
            return False
        
        protein = self._prot_service.annotate(tx_anno.protein_id)
        for feat in protein.protein_features:    
            if feat.info.name == self._feature and location.overlaps_with(feat.info.region):
                return True
        return False
    
    