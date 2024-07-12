from genophenocorr.model import Variant, VariantEffect, FeatureType
from genophenocorr.model.genome import Region
from ._api import VariantPredicate
from genophenocorr.preprocessing import ProteinMetadataService


class VariantEffectPredicate(VariantPredicate):
    """
    `VariantEffectPredicate` is a `VariantPredicate` that sets up testing
    for a specific `VariantEffect` on a given transcript ID. 
    
    Args:
        effect (VariantEffect): the variant effect to search for
        tx_id (str): the accession of the transcript of interest, e.g. `NM_123456.7`
    """
    
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
        
        tx_anno = variant.get_tx_anno_by_tx_id(self._tx_id)
        if tx_anno is None:
            return False
        for effect in tx_anno.variant_effects:
            if effect == self._effect:
                return True
        return False
    
    def __eq__(self, value: object) -> bool:
        if isinstance(value, VariantEffectPredicate):
            return self._effect == value._effect and self._tx_id == value._tx_id
        return False
    
    def __hash__(self) -> int:
        return hash((self._effect, self._tx_id))

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return f'VariantEffectPredicate(effect={self._effect}, tx_id={self._tx_id})'
    

class VariantKeyPredicate(VariantPredicate):
    """
    `VariantKeyPredicate` tests if the variant is described by the `key`.

    .. note::

      You can get the variant key by calling :class:`genophenocorr.model.VariantCoordinates.variant_key`.

    :param key: a `str` with the variant key.
    """
    
    def __init__(self, key: str) -> None:
        self._key = key
        
    def get_question(self) -> str:
        return f'variant has ID of {self._key}'

    def test(self, variant: Variant) -> bool:
        return self._key == variant.variant_coordinates.variant_key
    
    def __eq__(self, value: object) -> bool:
        if isinstance(value, VariantKeyPredicate):
            return self._key == value._key
        return False
    
    def __hash__(self) -> int:
        return hash((self._key,))

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return f'VariantKeyPredicate(key={self._key})'


class VariantGenePredicate(VariantPredicate):
    """
    `VariantGenePredicate` tests if the variant 
    is annotated to affect a gene denoted by a `symbol`.

    Args:
        symbol (str): the gene symbol, e.g. `FBN1`
    """
    
    def __init__(self, symbol: str) -> None:
        self._symbol = symbol

    def get_question(self) -> str:
        return f'variant affects gene {self._symbol}'

    def test(self, variant: Variant) -> bool:
        for tx in variant.tx_annotations:
            if tx.gene_id == self._symbol:
                return True
        return False
    
    def __eq__(self, value: object) -> bool:
        if isinstance(value, VariantGenePredicate):
            return self._symbol == value._symbol
        return False
    
    def __hash__(self) -> int:
        return hash((self._symbol,))

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return f'VariantGenePredicate(symbol={self._symbol})'


class VariantTranscriptPredicate(VariantPredicate):
    """
    `VariantTranscriptPredicate` tests if the variant 
    is annotated to affect a transcript with `tx_id` accession.

    Args:
        tx_id (str): the accessiono of the transcript of interest, e.g. `NM_123456.7`
    """
    
    def __init__(self, tx_id: str) -> None:
        self._tx_id = tx_id

    def get_question(self) -> str:
        return f'variant affects transcript {self._tx_id}'
        
    def test(self, variant: Variant) -> bool:
        for tx in variant.tx_annotations:
            if tx.transcript_id == self._tx_id:
                return True
        return False
    
    def __eq__(self, value: object) -> bool:
        if isinstance(value, VariantTranscriptPredicate):
            return self._tx_id == value._tx_id
        return False
    
    def __hash__(self) -> int:
        return hash((self._tx_id,))

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return f'VariantTranscriptPredicate(tx_id={self._tx_id})'
    
    
class VariantExonPredicate(VariantPredicate):
    """
    `VariantExonPredicate` tests if the variant affects 
    *n*-th exon of the transcript of interest.

    .. warning::

      We use 1-based numbering to number the exons, 
      not the usual 0-based numbering of the computer science.
      Therefore, the first exon of the transcript 
      has ``exon_number==1``, the second exon is ``2``, and so on ...

    .. warning::

      We do not check if the `exon_number` spans 
      beyond the number of exons of the given `transcript_id`!
      Therefore, ``exon_number==10,000`` will effectively 
      return :attr:`GenotypeBooleanPredicate.FALSE`
      for *all* patients!!! 😱
      Well, at least the patients of the *Homo sapiens sapiens* taxon...

    :param exon: a positive `int` of the target exon
    :param tx_id: the accession of the transcript of interest, e.g. `NM_123456.7`
    """
    
    def __init__(
        self,
        exon: int, 
        tx_id: str,
    ):
        self._exon = exon
        self._tx_id = tx_id
        
    def get_question(self) -> str:
        return f'variant affects exon {self._exon} on {self._tx_id}'

    def test(self, variant: Variant) -> bool:
        tx_anno = variant.get_tx_anno_by_tx_id(self._tx_id)
        if tx_anno is None:
            return False
        
        if tx_anno.overlapping_exons is None:
            return False

        return any(self._exon == exon for exon in tx_anno.overlapping_exons)
    
    def __eq__(self, value: object) -> bool:
        if isinstance(value, VariantExonPredicate):
            return self._tx_id == value._tx_id
        return False
    
    def __hash__(self) -> int:
        return hash((self._exon, self._tx_id,))

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return f'VariantExonPredicate(exon={self._exon}, tx_id={self._tx_id})'


class ProteinRegionPredicate(VariantPredicate):
    """
    `ProteinRegionPredicate` tests if the variant 
    overlaps with a given region of the protein.

    The predicate needs the start and end coordinate 
    for the protein region, given as a :class:`genophenocorr.model.genome.Region`.
    For instance, `Region(150, 175)`.

    :param region: a `Region` with the start and end coordinates
    :param tx_id: the accession of the transcript of interest
    """
    
    def __init__(
        self, 
        region: Region, 
        tx_id: str,
    ):
        self._region = region
        self._tx_id = tx_id
        
    def get_question(self) -> str:
        return f'variant affects aminoacid(s) between {self._region.start} and {self._region.end} on protein encoded by transcript {self._tx_id}'

    def test(self, variant: Variant) -> bool:
        tx_anno = variant.get_tx_anno_by_tx_id(self._tx_id)
        if tx_anno is None:
            return False
        location = tx_anno.protein_effect_location
        if location is None:
            return False
        return location.overlaps_with(self._region)
    
    def __eq__(self, value: object) -> bool:
        if isinstance(value, ProteinRegionPredicate):
            return self._region == value._region and self._tx_id == value._tx_id
        return False
    
    def __hash__(self) -> int:
        return hash((self._region, self._tx_id,))

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return f'ProteinRegionPredicate(region={self._region}, tx_id={self._tx_id})'

    
class ProteinFeatureTypePredicate(VariantPredicate):
    """
    `ProteinFeatureTypePredicate` checks if the variant overlaps with
    a protein :class:`genophenocorr.model.FeatureType` 
    in the protein encoded by a selected transcript.
    
    Args:
        feature_type (FeatureType): a member of the `FeatureType` enum
        tx_id (str): the accession of the selected transcript, e.g. `NM_123456.7`
        protein_metadata_service (ProteinMetadataService): a service for fetching protein information
    """
    
    def __init__(
        self, 
        feature_type: FeatureType, 
        tx_id: str, 
        protein_metadata_service: ProteinMetadataService,
    ):
        self._feature_type = feature_type
        self._tx_id = tx_id
        self._prot_service = protein_metadata_service
        
    def get_question(self) -> str:
        return f'variant affects {self._feature_type.name} feature type on the protein encoded by transcript {self._tx_id}'

    def test(self, variant: Variant) -> bool:
        tx_anno = variant.get_tx_anno_by_tx_id(self._tx_id)
        if tx_anno is None:
            return False
        location = tx_anno.protein_effect_location
        if location is None:
            return False
        
        protein_id = tx_anno.protein_id
        if protein_id is None:
            return False
        
        protein_meta = self._prot_service.annotate(protein_id)
        for feat in protein_meta.protein_features:    
            if feat.feature_type == self._feature_type and location.overlaps_with(feat.info.region):
                return True
        return False

    def __eq__(self, value: object) -> bool:
        if isinstance(value, ProteinFeatureTypePredicate):
            # We do not care about `self._prot_service`
            return self._feature_type == value._feature_type and self._tx_id == value._tx_id
        return False
    
    def __hash__(self) -> int:
        # We do not care about `self._prot_service`
        return hash((self._feature_type, self._tx_id,))

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return f'ProteinFeatureTypePredicate(feature_type={self._feature_type}, tx_id={self._tx_id}, protein_metadata_service={self._prot_service})'


class ProteinFeaturePredicate(VariantPredicate):
    """
    `ProteinFeaturePredicate` checks if the variant overlaps with
    a feature type in the protein encoded by a transcript of interest.
    
    Args:
        feature_type (FeatureType): a member of the `FeatureType` enum
        tx_id (str): the accession of the transcript of interest, e.g. `NM_123456.7`
        protein_metadata_service (ProteinMetadataService): a service for fetching protein information
    """
    
    def __init__(self, feature_name: str, tx_id: str, protein_metadata_service: ProteinMetadataService) -> None:
        self._feature_name = feature_name
        self._tx_id = tx_id
        self._prot_service = protein_metadata_service

    def get_question(self) -> str:
        return f'Variant that affects {self._feature_name} feature on the protein encoded by transcript {self._tx_id}'
        
    def test(self, variant: Variant) -> bool:
        tx_anno = variant.get_tx_anno_by_tx_id(self._tx_id)
        if tx_anno is None:
            return False
        
        location = tx_anno.protein_effect_location
        if location is None:
            return False
        
        protein_id = tx_anno.protein_id
        if protein_id is None:
            return False

        protein = self._prot_service.annotate(protein_id)
        for feat in protein.protein_features:    
            if feat.info.name == self._feature_name and location.overlaps_with(feat.info.region):
                return True
        return False
    
    def __eq__(self, value: object) -> bool:
        if isinstance(value, ProteinFeaturePredicate):
            # We do not care about `self._prot_service`
            return self._feature_name == value._feature_name and self._tx_id == value._tx_id
        return False
    
    def __hash__(self) -> int:
        # We do not care about `self._prot_service`
        return hash((self._feature_name, self._tx_id,))

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return f'ProteinFeaturePredicate(feature_name={self._feature_name}, tx_id={self._tx_id}, protein_metadata_service={self._prot_service})'
