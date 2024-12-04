import typing
import hpotk

import operator

from gpsea.model import Variant, VariantEffect, VariantClass, ProteinMetadata, FeatureType
from gpsea.model.genome import Region

from ._api import VariantPredicate


class AlwaysTrueVariantPredicate(VariantPredicate):
    """
    `AlwaysTrueVariantPredicate` returns `True` for any variant.
    """

    @staticmethod
    def get_instance() -> "AlwaysTrueVariantPredicate":
        return ALWAYS_TRUE

    @property
    def name(self) -> str:
        return "Always True"

    @property
    def description(self) -> str:
        return "true"

    @property
    def variable_name(self) -> str:
        return "N/A"

    def test(self, variant: Variant) -> bool:
        return True

    def __eq__(self, value: object) -> bool:
        return isinstance(value, AlwaysTrueVariantPredicate)
    
    def __hash__(self) -> int:
        return 17

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return 'AlwaysTrueVariantPredicate()'


ALWAYS_TRUE = AlwaysTrueVariantPredicate()


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
    
    @property
    def name(self) -> str:
        return "Variant Effect"

    @property
    def description(self) -> str:
        return f"{self._effect.name} on {self._tx_id}"

    @property
    def variable_name(self) -> str:
        return "variant effect"

    def test(self, variant: Variant) -> bool:
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

      You can get the variant key by calling :class:`~gpsea.model.VariantCoordinates.variant_key`.

    :param key: a `str` with the variant key.
    """
    
    def __init__(self, key: str) -> None:
        self._key = key

    @property
    def name(self) -> str:
        return "Variant Key Predicate"

    @property
    def description(self) -> str:
        return f"variant key is {self._key}"

    @property
    def variable_name(self) -> str:
        return "variant key"

    def test(self, variant: Variant) -> bool:
        return self._key == variant.variant_info.variant_key
    
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

    @property
    def name(self) -> str:
        return "Gene Predicate"

    @property
    def description(self) -> str:
        return f"affects {self._symbol}"

    @property
    def variable_name(self) -> str:
        return "gene"

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

    @property
    def name(self) -> str:
        return "Transcript Predicate"

    @property
    def description(self) -> str:
        return f"affects {self._tx_id}"

    @property
    def variable_name(self) -> str:
        return "transcript"
        
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
      return `False` for *all* variants!!! 😱
      Well, at least the genomic variants of the *Homo sapiens sapiens* taxon...

    :param exon: a positive `int` of the target exon
    :param tx_id: the accession of the transcript of interest, e.g. `NM_123456.7`
    """
    
    def __init__(
        self,
        exon: int,
        tx_id: str,
    ):
        assert isinstance(exon, int) and exon > 0, '`exon` must be a positive `int`'
        self._exon = exon
        assert isinstance(tx_id, str)
        self._tx_id = tx_id
    
    @property
    def name(self) -> str:
        return "Exon Predicate"

    @property
    def description(self) -> str:
        return f"overlaps with exon {self._exon} of {self._tx_id}"

    @property
    def variable_name(self) -> str:
        return "exon"

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


class IsLargeImpreciseStructuralVariantPredicate(VariantPredicate):
    """
    `IsLargeImpreciseStructuralVariantPredicate` tests if the variant is a large imprecise SV
    where the exact breakpoint coordinates are not available.
    """

    @property
    def name(self) -> str:
        return "Large Imprecise SV"

    @property
    def description(self) -> str:
        return "is large imprecise structural variant"

    @property
    def variable_name(self) -> str:
        return "variant"

    def test(self, variant: Variant) -> bool:
        return variant.variant_info.has_sv_info()

    def __eq__(self, value: object) -> bool:
        return isinstance(value, IsLargeImpreciseStructuralVariantPredicate)
    
    def __hash__(self) -> int:
        return hash(())

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return 'IsLargeImpreciseStructuralVariantPredicate'


class VariantClassPredicate(VariantPredicate):
    """
    `VariantClassPredicate` tests if the variant class matches the query (e.g. `DEL`).
    """

    def __init__(
        self,
        query: VariantClass
    ):
        assert isinstance(query, VariantClass), 'query must be `VariantClass`'
        self._query = query
    
    @property
    def name(self) -> str:
        return "Variant Class"

    @property
    def description(self) -> str:
        return f"variant class is {self._query.name}"

    @property
    def variable_name(self) -> str:
        return "variant class"

    def test(self, variant: Variant) -> bool:
        """
        We are testing for a structural variant that is a deletion.
        Deletions can be represented in two ways in our data.
        1. Imprecisely, using the sequence ontology term for deletion (see code)
        2. Using the VCF-syntax "DEL"
        """
        return variant.variant_info.variant_class == self._query

    def __eq__(self, value: object) -> bool:
        return isinstance(value, VariantClassPredicate) and self._query == value._query
    
    def __hash__(self) -> int:
        return hash((self._query,))

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return f'VariantClassPredicate(query={self._query})'
    

class StructuralTypePredicate(VariantPredicate):
    """
    `StructuralTypePredicate` tests if the variant has a certain structural type,
    e.g. `chromosomal_deletion` (`SO:1000029 <http://purl.obolibrary.org/obo/SO_1000029>`_).
    """

    @staticmethod
    def from_curie(curie: typing.Union[str, hpotk.TermId]):
        if isinstance(curie, str):
            query = hpotk.TermId.from_curie(curie)
        elif isinstance(curie, hpotk.TermId):
            query = curie
        else:
            raise ValueError(f'curie `{curie}` must be a `str` or `TermId` but was f{type(curie)}')
        
        return StructuralTypePredicate(query)

    def __init__(
        self,
        query: hpotk.TermId,
    ):
        assert isinstance(query, hpotk.TermId), 'query must be a `TermId`'
        self._query = query

    @property
    def name(self) -> str:
        return "Structural Type"

    @property
    def description(self) -> str:
        return f"structural type is {self._query.value}"

    @property
    def variable_name(self) -> str:
        return "structural type"

    def test(self, variant: Variant) -> bool:
        sv_info = variant.variant_info.sv_info
        if sv_info is not None:
            return sv_info.structural_type == self._query
        return False

    def __eq__(self, value: object) -> bool:
        return isinstance(value, StructuralTypePredicate) and self._query == value._query
    
    def __hash__(self) -> int:
        return hash((self._query,))

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return f'StructuralTypePredicate(query={self._query})'


def _decode_operator(op: str) -> typing.Callable[[int, int], bool]:
    if op == '<':
        return operator.lt
    elif op == '<=':
        return operator.le
    elif op == '==':
        return operator.eq
    elif op == '!=':
        return operator.ne
    elif op == '>=':
        return operator.ge
    elif op == '>':
        return operator.gt
    else:
        raise ValueError(f'Unsupported operator {op}')


class ChangeLengthPredicate(VariantPredicate):
    """
    `ChangeLengthPredicate` tests if the variant's change length is above, below, or equal to a threshold.
    """
    
    def __init__(
        self,
        operator: typing.Literal['<', '<=', '==', '!=', '>=', '>'],
        threshold: int,
    ):
        self._operator_str = operator
        self._operator = _decode_operator(operator)
        
        assert isinstance(threshold, int), 'threshold must be an `int`'
        self._threshold = threshold

    @property
    def name(self) -> str:
        return "Change Length"

    @property
    def description(self) -> str:
        return f"change length {self._operator_str} {self._threshold}"

    @property
    def variable_name(self) -> str:
        return "change length"

    def test(self, variant: Variant) -> bool:
        vc = variant.variant_info.variant_coordinates
        if vc is not None:
            return self._operator(vc.change_length, self._threshold)
        return False

    def __eq__(self, value: object) -> bool:
        return isinstance(value, ChangeLengthPredicate) \
            and self._operator_str == value._operator_str \
            and self._threshold == value._threshold
    
    def __hash__(self) -> int:
        return hash((self._operator_str, self._threshold,))

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return f'ChangeLengthPredicate(operator=\'{self._operator_str}\', threshold={self._threshold})'


class RefAlleleLengthPredicate(VariantPredicate):
    """
    `RefAlleleLengthPredicate` tests if the length of the variant's reference allele
    is greater than, equal, or less than certain value.
    """

    def __init__(
        self,
        operator: typing.Literal['<', '<=', '==', '!=', '>=', '>'],
        length: int,
    ):
        self._operator_str = operator
        self._operator = _decode_operator(operator)
        
        assert isinstance(length, int), 'length must be an `int`'
        assert length >= 0, 'length must be non-negative'
        self._length = length

    @property
    def name(self) -> str:
        return "Reference Allele Length"

    @property
    def description(self) -> str:
        return f"reference allele length {self._operator_str} {self._length}"

    @property
    def variable_name(self) -> str:
        return "length of the reference allele"

    def test(self, variant: Variant) -> bool:
        vc = variant.variant_info.variant_coordinates
        if vc is not None:
            return self._operator(len(vc), self._length)
        return False

    def __eq__(self, value: object) -> bool:
        return isinstance(value, RefAlleleLengthPredicate) \
            and self._operator_str == value._operator_str \
            and self._length == value._length
    
    def __hash__(self) -> int:
        return hash((self._operator_str, self._length,))

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return f'RefAlleleLengthPredicate(operator=\'{self._operator_str}\', length={self._length})'


class ProteinRegionPredicate(VariantPredicate):
    """
    `ProteinRegionPredicate` tests if the variant
    overlaps with a given region of the protein.

    The predicate needs the start and end coordinate
    for the protein region, given as a :class:`~gpsea.model.genome.Region`.
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
        
    @property
    def name(self) -> str:
        return "Protein Region"

    @property
    def description(self) -> str:
        return (
            "overlaps with "
            f"[{self._region.start + 1},{self._region.end}] region "  # Reporting 1-based coordinates
            f"of the protein encoded by {self._tx_id}"
        )

    @property
    def variable_name(self) -> str:
        return "overlap with aminoacid region"

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
    a protein :class:`~gpsea.model.FeatureType`
    in the protein encoded by a selected transcript.
    
    Args:
        feature_type (FeatureType): a member of the `FeatureType` enum
        tx_id (str): the accession of the selected transcript, e.g. `NM_123456.7`
        protein_metadata_service (ProteinMetadataService): a service for fetching protein information
    """
    
    def __init__(
        self,
        feature_type: typing.Union[FeatureType, str],
        protein_metadata: ProteinMetadata,
    ):
        if isinstance(feature_type, str):
            self._feature_type = feature_type
        elif isinstance(feature_type, FeatureType):
            FeatureType.deprecation_warning()
            self._feature_type = feature_type.name
        else:
            raise ValueError(f'`feature_type` must be either `FeatureType` or `str` but was {type(feature_type)}')
        
        assert isinstance(protein_metadata, ProteinMetadata)
        self._protein_metadata = protein_metadata
        
    @property
    def name(self) -> str:
        return "Protein Feature Type"

    @property
    def description(self) -> str:
        return f"overlaps with a protein feature {self._feature_type}"

    @property
    def variable_name(self) -> str:
        return f"overlap with {self._feature_type} feature"

    def test(self, variant: Variant) -> bool:
        for tx_ann in variant.tx_annotations:
            if tx_ann.protein_id is not None and tx_ann.protein_id == self._protein_metadata.protein_id:
                location = tx_ann.protein_effect_location
                if location is not None:
                    return any(
                        self._feature_type == feature.feature_type and location.overlaps_with(feature.info.region)
                        for feature in self._protein_metadata.protein_features
                    )
                
        return False

    def __eq__(self, value: object) -> bool:
        return isinstance(value, ProteinFeatureTypePredicate) \
            and self._feature_type == value._feature_type \
            and self._protein_metadata == value._protein_metadata
    
    def __hash__(self) -> int:
        return hash((self._feature_type, self._protein_metadata,))

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return 'ProteinFeatureTypePredicate(' \
            f'feature_type={self._feature_type}, ' \
            f'protein_metadata={self._protein_metadata})'


class ProteinFeaturePredicate(VariantPredicate):
    """
    `ProteinFeaturePredicate` checks if the variant overlaps with
    a feature with a specific id.
    
    Args:
        feature_type (FeatureType): a member of the `FeatureType` enum
        tx_id (str): the accession of the transcript of interest, e.g. `NM_123456.7`
        protein_metadata_service (ProteinMetadataService): a service for fetching protein information
    """
    
    def __init__(
        self,
        feature_id: str,
        protein_metadata: ProteinMetadata,
    ):
        assert isinstance(feature_id, str)
        self._feature_id = feature_id
        
        assert isinstance(protein_metadata, ProteinMetadata)
        self._protein_metadata = protein_metadata

    @property
    def name(self) -> str:
        return "Protein Feature"

    @property
    def description(self) -> str:
        return f"overlaps with {self._feature_id}"

    @property
    def variable_name(self) -> str:
        return "overlap with a protein feature"
        
    def test(self, variant: Variant) -> bool:
        for tx_ann in variant.tx_annotations:
            if tx_ann.protein_id is not None and tx_ann.protein_id == self._protein_metadata.protein_id:
                location = tx_ann.protein_effect_location
                if location is not None:
                    return any(
                        self._feature_id == feature.info.name and location.overlaps_with(feature.info.region)
                        for feature in self._protein_metadata.protein_features
                    )
        
        return False
    
    def __eq__(self, value: object) -> bool:
        return isinstance(value, ProteinFeaturePredicate) \
            and self._feature_id == value._feature_id \
            and self._protein_metadata == value._protein_metadata
    
    def __hash__(self) -> int:
        # We do not care about `self._prot_service`
        return hash((self._feature_id, self._protein_metadata,))

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return 'ProteinFeaturePredicate(' \
            f'feature_id={self._feature_id}, ' \
            f'protein_metadata={self._protein_metadata})'
