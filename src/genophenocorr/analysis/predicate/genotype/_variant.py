import abc
import typing

from genophenocorr.model import VariantEffect, Variant, FeatureType
from genophenocorr.preprocessing import ProteinMetadataService
from ._api import VariantPredicate


class LogicalVariantPredicate(VariantPredicate, metaclass=abc.ABCMeta):
    # NOT PART OF THE PUBLIC API

    def __init__(
            self,
            predicates: typing.Iterable[VariantPredicate],
    ):
        self._predicates = tuple(predicates)


class AnyVariantPredicate(LogicalVariantPredicate):
    # NOT PART OF THE PUBLIC API

    def test(self, variant: Variant) -> bool:
        return any(predicate.test(variant) for predicate in self._predicates)


class AllVariantPredicate(LogicalVariantPredicate):
    # NOT PART OF THE PUBLIC API

    def test(self, variant: Variant) -> bool:
        return all(predicate.test(variant) for predicate in self._predicates)


class VariantPredicates:
    """
    `VariantPredicates` is a static utility class to provide the variant predicates
    that are relatively simple to configure.
    """

    @staticmethod
    def variant_effect(effect: VariantEffect) -> VariantPredicate:
        """
        Prepare a :class:`VariantPredicate` to test if the functional annotation predicts the variant to lead to
        a certain variant effect.
        Args:
            effect: the target :class:`VariantEffect`

        Returns:
            VariantPredicate: a predicate for testing
        """
        raise NotImplementedError('Not yet implemented')

    @staticmethod
    def variant_key(key: str) -> VariantPredicate:
        """
        Prepare a :class:`VariantPredicate` that tests if the variant matches the provided `key`.
        Args:
            key: a `str` with the variant key (e.g. `X_12345_12345_C_G` or `22_10001_20000_INV`)

        Returns:
            VariantPredicate: a predicate
        """
        raise NotImplementedError('Not yet implemented')

    @staticmethod
    def gene(symbol: str) -> VariantPredicate:
        """
        Prepare a :class:`VariantPredicate` that tests if the variant affects a given gene.
        Args:
            symbol: a `str` with the gene symbol (e.g. `FBN1`).

        Returns:
            VariantPredicate: a predicate
        """

    @staticmethod
    def transcript(tx_id: str) -> VariantPredicate:
        """
        Prepare a :class:`VariantPredicate` that tests if the variant affects a transcript.
        Args:
            tx_id: a `str` with the transcript accession.

        Returns:
            VariantPredicate: a predicate
        """

    @staticmethod
    def exon(exon: int) -> VariantPredicate:
        """
        Prepare a :class:`VariantPredicate` that tests if the variant overlaps with an exon of a specific transcript.
        Args:
            exon: a non-negative `int` with the index of the target exon (e.g. `0` for the 1st exon, `1` for the 2nd, ...)

        Returns:
            VariantPredicate: a predicate
        """
        raise NotImplementedError('Not yet implemented')

    @staticmethod
    def all_match(predicates: typing.Iterable[VariantPredicate]) -> VariantPredicate:
        """
        Prepare a :class:`VariantPredicate` that tests if the variant passes testing
        with *all* of the provided `predicates`.
        Args:
            predicates: an iterable with the :class:`VariantPredicate` instances to test.

        Returns:
            VariantPredicate: a predicate
        """
        return AllVariantPredicate(predicates)

    @staticmethod
    def any_match(predicates: typing.Iterable[VariantPredicate]) -> VariantPredicate:
        """
        Prepare a :class:`VariantPredicate` that tests if the variant passes testing
        with *any* of the provided `predicates`.
        Args:
            predicates: an iterable with the :class:`VariantPredicate` instances to test.

        Returns:
            VariantPredicate: a predicate
        """
        return AnyVariantPredicate(predicates)


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
    ) -> VariantPredicate:
        """
        Prepare a :class:`VariantPredicate` that tests if the variant affects a protein feature type.
        Args:
            feature_type: the target protein :class:`FeatureType` (e.g. ``)

        Returns:
            VariantPredicate: a predicate
        """
        raise NotImplementedError('Not yet implemented')

    def protein_feature(
            self,
            feature_id: str,
    ) -> VariantPredicate:
        """
        Prepare a :class:`VariantPredicate` that tests if the variant affects a protein feature type.
        Args:
            feature_id: the id of the target protein feature (e.g. `ANK 1`)

        Returns:
            VariantPredicate: a predicate
        """
        raise NotImplementedError('Not yet implemented')
