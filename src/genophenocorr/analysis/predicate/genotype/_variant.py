import typing

from genophenocorr.model import VariantEffect, FeatureType
from genophenocorr.model.genome import Region
from genophenocorr.preprocessing import ProteinMetadataService
from ._api import VariantPredicate
from ._predicates import *


# We do not need more than just one instance of these predicates.
IS_BND = VariantClassPredicate(VariantClass.BND)
IS_LARGE_IMPRECISE_SV = IsLargeImpreciseStructuralVariantPredicate()


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
    def region(region: Region, tx_id: str) -> VariantPredicate:
        """
        Prepare a :class:`VariantPredicate` that tests if the variant overlaps with a region on a protein of a specific transcript.

        Args:
            region: a :class:`Region` that gives the start and end coordinate of the region of interest on a protein strand.

        Returns:
            VariantPredicate: a predicate
        """
        return ProteinRegionPredicate(region, tx_id)

    @staticmethod
    def is_large_imprecise_sv() -> VariantPredicate:
        """
        Prepare a :class:`VariantPredicate` for testing if the variant is a large structural variant (SV)
        without exact breakpoint coordinates.
        """
        return IS_LARGE_IMPRECISE_SV

    @staticmethod
    def is_structural_variant(
        threshold: int = 50,
    ) -> VariantPredicate:
        """
        Prepare a :class:`VariantPredicate` for testing if the variant is a structural variant (SV).

        SVs are usually defined as variant affecting more than a certain number of base pairs.
        The thresholds vary in the literature, but here we use 50bp as a default.

        Any variant that affects at least `threshold` base pairs is considered an SV. 
        Large SVs with unknown breakpoint coordinates or translocations (:class:`VariantClass.BND`) 
        are always considered as an SV.

        Args:
            threshold: a non-negative `int` with the number of base pairs that must be affected
        """
        assert threshold >= 0, '`threshold` must be non-negative!'
        return VariantPredicates.change_length('<=', -threshold) | VariantPredicates.change_length('>=', threshold) | VariantPredicates.is_large_imprecise_sv() | IS_BND

    @staticmethod
    def structural_type(
        curie: typing.Union[str, hpotk.TermId],
    ) -> VariantPredicate:
        """
        Prepare a :class:`VariantPredicate` for testing if the variant has a certain structural type.

        We recommend using a descendant of `structural_variant` (`SO:0001537 <https://purl.obolibrary.org/obo/SO_0001537>`_)
        as the structural type.

        **Example**

        Make a predicate for testing if the variant is a chromosomal deletion (`SO:1000029`):

        >>> from genophenocorr.analysis.predicate.genotype import VariantPredicates
        >>> predicate = VariantPredicates.structural_type('SO:1000029')
        >>> predicate
        StructuralTypePredicate(query=SO:1000029)

        Args:
            curie: compact uniform resource identifier (CURIE) with the structural type to test.

        Returns:
            VariantPredicate: a predicate
        """
        return StructuralTypePredicate.from_curie(curie)

    @staticmethod
    def variant_class(
        variant_class: VariantClass,
    ) -> VariantPredicate:
        """
        Prepare a :class:`VariantPredicate` for testing if the variant is of a certain :class:`VariantClass`.

        Args:
            variant_class: the variant class to test.

        Returns:
            VariantPredicate: a predicate
        """
        return VariantClassPredicate(
            query=variant_class,
        )

    @staticmethod
    def change_length(
        operator: typing.Literal["<", "<=", "==", "!=", ">=", ">"],
        threshold: int,
    ) -> VariantPredicate:
        """
        Prepare a :class:`VariantPredicate` for testing if the variant's change length
        is above, below, or (not) equal to certain `threshold`.

        .. seealso::

            See :meth:`genophenocorr.model.VariantCoordinates.change_length` for more info on change length.

        **Example**

        Make a predicate for testing if the change length is less than or equal to `-10`, 
        e.g. to test if a variant is a *deletion* leading to removal of at least 10 base pairs:

        >>> from genophenocorr.analysis.predicate.genotype import VariantPredicates
        >>> predicate = VariantPredicates.change_length('<=', -10)
        >>> predicate
        ChangeLengthPredicate(operator='<=', threshold=-10)

        Args:
            operator: a `str` with the desired test. Must be one of ``{ '<', '<=', '==', '!=', '>=', '>' }``.
        """
        return ChangeLengthPredicate(operator, threshold)

    @staticmethod
    def is_structural_deletion(
        threshold: int = -50,
    ) -> VariantPredicate:
        """
        Prepare a :class:`VariantPredicate` for testing if the variant 
        is a `chromosomal deletion <https://purl.obolibrary.org/obo/SO_1000029>`_ or a structural variant deletion 
        that leads to removal of at least *n* base pairs (50bp by default).
        
        .. note::

            The predicate uses See :meth:`genophenocorr.model.VariantCoordinates.change_length`
            to determine if the length of the variant is above or below `threshold`.
            
            **IMPORTANT**: the change lengths of deletions are *negative*, since the alternate allele is shorter than the reference allele.

        **Example**

        Prepare a predicate for testing if the variant is a chromosomal destructural deletion with more than 

        Args:
            threshold: an `int` with the change length threshold to determine if a variant is "structural" (-50 bp by default).

        Returns:
            VariantPredicate: a predicate
        """
        chromosomal_deletion = "SO:1000029"
        return VariantPredicates.structural_type(chromosomal_deletion) | (
            VariantPredicates.variant_class(VariantClass.DEL)
            & VariantPredicates.change_length("<=", threshold)
        )


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
        return ProteinFeatureTypePredicate(
            feature_type, tx_id, self._protein_metadata_service
        )

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
        return ProteinFeaturePredicate(
            feature_id, tx_id, self._protein_metadata_service
        )
