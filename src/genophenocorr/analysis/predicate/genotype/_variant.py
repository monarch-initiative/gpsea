from genophenocorr.model import VariantEffect, FeatureType
from genophenocorr.model.genome import Region
from genophenocorr.preprocessing import ProteinMetadataService
from ._api import VariantPredicate
from ._predicates import *

class VariantPredicates:
    """
    `VariantPredicates` is a static utility class to provide the variant predicates
    that are relatively simple to configure.
    """

    @staticmethod
    def variant_effect(
        effect: VariantEffect,
        tx_id: str,
    ) -> VariantPredicate:
        """
        Prepare a :class:`VariantPredicate` to test if the functional annotation predicts the variant to lead to
        a certain variant effect.
        Args:
            effect: the target :class:`VariantEffect`
            tx_id: a `str` with the accession ID of the target transcript (e.g. `NM_123.4`)

        Returns:
            VariantPredicate: a predicate for testing
        """
        return VariantEffectPredicate(effect, tx_id)

    @staticmethod
    def variant_key(key: str) -> VariantPredicate:
        """
        Prepare a :class:`VariantPredicate` that tests if the variant matches the provided `key`.
        Args:
            key: a `str` with the variant key (e.g. `X_12345_12345_C_G` or `22_10001_20000_INV`)

        Returns:
            VariantPredicate: a predicate
        """
        return VariantKeyPredicate(key)

    @staticmethod
    def gene(symbol: str) -> VariantPredicate:
        """
        Prepare a :class:`VariantPredicate` that tests if the variant affects a given gene.
        Args:
            symbol: a `str` with the gene symbol (e.g. `FBN1`).

        Returns:
            VariantPredicate: a predicate
        """
        return VariantGenePredicate(symbol)

    @staticmethod
    def transcript(tx_id: str) -> VariantPredicate:
        """
        Prepare a :class:`VariantPredicate` that tests if the variant affects a transcript.
        Args:
            tx_id: a `str` with the accession ID of the target transcript (e.g. `NM_123.4`)

        Returns:
            VariantPredicate: a predicate
        """
        return VariantTranscriptPredicate(tx_id)

    @staticmethod
    def exon(
        exon: int,
        tx_id: str,
    ) -> VariantPredicate:
        """
        Prepare a :class:`VariantPredicate` that tests if the variant overlaps with an exon of a specific transcript.
        Args:
            exon: a non-negative `int` with the index of the target exon (e.g. `0` for the 1st exon, `1` for the 2nd, ...)
            tx_id: a `str` with the accession ID of the target transcript (e.g. `NM_123.4`)

        Returns:
            VariantPredicate: a predicate
        """
        return VariantExonPredicate(exon, tx_id)
    
    @staticmethod
    def region(
        region: Region, 
        tx_id: str
    ) -> VariantPredicate:
        """
        Prepare a :class:`VariantPredicate` that tests if the variant overlaps with a region on a protein of a specific transcript.
        Args:
            region: a :class:`Region` that gives the start and end coordinate of the region of interest on a protein strand.

        Returns:
            VariantPredicate: a predicate
        """
        return ProteinRegionPredicate(region, tx_id)
    
    @staticmethod
    def is_structural_variant() -> VariantPredicate:
        """
        Prepare a :class:`VariantPredicate` for testing if the variant is a structural variant.

        Check :meth:`genophenocorr.model.VariantInfo.is_structural` for the 

        Returns:
            VariantPredicate: a predicate
        """
        return IS_STRUCTURAL


class ProteinPredicates:
    """
    `ProteinPredicates` prepares variant predicates that need to consult :class:`ProteinMetadataService`
    to categorize a :class:`Variant`.
    """

    def __init__(
            self,
            protein_metadata_service: ProteinMetadataService,
    ):
        self._protein_metadata_service = protein_metadata_service

    def protein_feature_type(
            self,
            feature_type: FeatureType,
            tx_id: str
    ) -> VariantPredicate:
        """
        Prepare a :class:`VariantPredicate` that tests if the variant affects a protein feature type.
        Args:
            feature_type: the target protein :class:`FeatureType` (e.g. ``)

        Returns:
            VariantPredicate: a predicate
        """
        return ProteinFeatureTypePredicate(feature_type, tx_id, self._protein_metadata_service)

    def protein_feature(
            self,
            feature_id: str,
            tx_id: str
    ) -> VariantPredicate:
        """
        Prepare a :class:`VariantPredicate` that tests if the variant affects a protein feature type.
        Args:
            feature_id: the id of the target protein feature (e.g. `ANK 1`)

        Returns:
            VariantPredicate: a predicate
        """
        return ProteinFeaturePredicate(feature_id, tx_id, self._protein_metadata_service)
