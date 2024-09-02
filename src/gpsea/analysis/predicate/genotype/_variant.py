import typing

from gpsea.model import VariantEffect, FeatureType
from gpsea.model.genome import Region
from gpsea.preprocessing import ProteinMetadataService
from ._api import VariantPredicate, AllVariantPredicate, AnyVariantPredicate
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
    def all(predicates: typing.Iterable[VariantPredicate]) -> VariantPredicate:
        """
        Prepare a :class:`VariantPredicate` that returns `True` if ALL `predicates` evaluate to `True`.

        This is useful for building compound predicates programmatically.

        **Example**

        Build a predicate to test if variant has a functional annotation to genes *SURF1* and *SURF2*:

        >>> from gpsea.analysis.predicate.genotype import VariantPredicates

        >>> genes = ('SURF1', 'SURF2',)
        >>> predicate = VariantPredicates.all(VariantPredicates.gene(g) for g in genes)
        >>> predicate.get_question()
        '(impacts SURF1 AND impacts SURF2)'

        Args:
            predicates: an iterable of predicates to test
        """
        predicates = tuple(predicates)
        if len(predicates) == 1:
            # No need to wrap one predicate into a logical predicate.
            return predicates[0]
        else:
            return AllVariantPredicate(*predicates)

    @staticmethod
    def any(predicates: typing.Iterable[VariantPredicate]) -> VariantPredicate:
        """
        Prepare a :class:`VariantPredicate` that returns `True` if ANY of the `predicates` evaluates to `True`.

        This can be useful for building compound predicates programmatically.

        **Example**

        Build a predicate to test if variant leads to a missense
        or nonsense change on a fictional transcript `NM_123456.7`:

        >>> from gpsea.model import VariantEffect
        >>> from gpsea.analysis.predicate.genotype import VariantPredicates

        >>> tx_id = 'NM_123456.7'
        >>> effects = (VariantEffect.MISSENSE_VARIANT, VariantEffect.STOP_GAINED,)
        >>> predicate = VariantPredicates.any(VariantPredicates.variant_effect(e, tx_id) for e in effects)
        >>> predicate.get_question()
        '(MISSENSE_VARIANT on NM_123456.7 OR STOP_GAINED on NM_123456.7)'

        Args:
            predicates: an iterable of predicates to test
        """
        predicates = tuple(predicates)
        if len(predicates) == 1:
            # No need to wrap one predicate into a logical predicate.
            return predicates[0]
        else:
            return AnyVariantPredicate(*predicates)

    @staticmethod
    def variant_effect(
        effect: VariantEffect,
        tx_id: str,
    ) -> VariantPredicate:
        """
        Prepare a :class:`VariantPredicate` to test if the functional annotation predicts the variant to lead to
        a certain variant effect.

        **Example**

        Make a predicate for testing if the variant leads to a missense change on transcript `NM_123.4`:

        >>> from gpsea.model import VariantEffect
        >>> from gpsea.analysis.predicate.genotype import VariantPredicates
        >>> predicate = VariantPredicates.variant_effect(VariantEffect.MISSENSE_VARIANT, tx_id='NM_123.4')
        >>> predicate.get_question()
        'MISSENSE_VARIANT on NM_123.4'

        Args:
            effect: the target :class:`VariantEffect`
            tx_id: a `str` with the accession ID of the target transcript (e.g. `NM_123.4`)
        """
        return VariantEffectPredicate(effect, tx_id)

    @staticmethod
    def variant_key(key: str) -> VariantPredicate:
        """
        Prepare a :class:`VariantPredicate` that tests if the variant matches the provided `key`.

        Args:
            key: a `str` with the variant key (e.g. `X_12345_12345_C_G` or `22_10001_20000_INV`)
        """
        return VariantKeyPredicate(key)

    @staticmethod
    def gene(symbol: str) -> VariantPredicate:
        """
        Prepare a :class:`VariantPredicate` that tests if the variant affects a given gene.

        Args:
            symbol: a `str` with the gene symbol (e.g. ``'FBN1'``).
        """
        return VariantGenePredicate(symbol)

    @staticmethod
    def transcript(tx_id: str) -> VariantPredicate:
        """
        Prepare a :class:`VariantPredicate` that tests if the variant affects a transcript.

        Args:
            tx_id: a `str` with the accession ID of the target transcript (e.g. `NM_123.4`)
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
            exon: a non-negative `int` with the index of the target exon
                (e.g. `0` for the 1st exon, `1` for the 2nd, ...)
            tx_id: a `str` with the accession ID of the target transcript (e.g. `NM_123.4`)
        """
        return VariantExonPredicate(exon, tx_id)

    @staticmethod
    def region(region: Region, tx_id: str) -> VariantPredicate:
        """
        Prepare a :class:`VariantPredicate` that tests if the variant
        overlaps with a region on a protein of a specific transcript.

        Args:
            region: a :class:`Region` that gives the start and end coordinate
                of the region of interest on a protein strand.
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
        assert threshold >= 0, "`threshold` must be non-negative!"
        return (
            VariantPredicates.change_length("<=", -threshold)
            | VariantPredicates.change_length(">=", threshold)
            | VariantPredicates.is_large_imprecise_sv()
            | IS_BND
        )

    @staticmethod
    def structural_type(
        curie: typing.Union[str, hpotk.TermId],
    ) -> VariantPredicate:
        """
        Prepare a :class:`VariantPredicate` for testing if the variant has a certain structural type.

        We recommend using a descendant of `structural_variant`
        (`SO:0001537 <https://purl.obolibrary.org/obo/SO_0001537>`_)
        as the structural type.

        **Example**

        Make a predicate for testing if the variant is a chromosomal deletion (`SO:1000029`):

        >>> from gpsea.analysis.predicate.genotype import VariantPredicates
        >>> predicate = VariantPredicates.structural_type('SO:1000029')
        >>> predicate.get_question()
        'structural type is SO:1000029'

        Args:
            curie: compact uniform resource identifier (CURIE) with the structural type to test.
        """
        return StructuralTypePredicate.from_curie(curie)

    @staticmethod
    def variant_class(
        variant_class: VariantClass,
    ) -> VariantPredicate:
        """
        Prepare a :class:`VariantPredicate` for testing if the variant is of a certain :class:`VariantClass`.

        **Example**

        Make a predicate to test if the variant is a deletion:

        >>> from gpsea.model import VariantClass
        >>> from gpsea.analysis.predicate.genotype import VariantPredicates
        >>> predicate = VariantPredicates.variant_class(VariantClass.DEL)
        >>> predicate.get_question()
        'variant class is DEL'

        Args:
            variant_class: the variant class to test.
        """
        return VariantClassPredicate(
            query=variant_class,
        )

    @staticmethod
    def ref_length(
        operator: typing.Literal["<", "<=", "==", "!=", ">=", ">"],
        length: int,
    ) -> VariantPredicate:
        """
        Prepare a :class:`VariantPredicate` for testing if the reference (REF) allele
        of variant is above, below, or (not) equal to certain `length`.

        .. seealso::

            See :ref:`length-of-the-reference-allele` for more info.

        **Example**

        Prepare a predicate that tests that the REF allele includes more than 5 base pairs:

        >>> from gpsea.analysis.predicate.genotype import VariantPredicates
        >>> predicate = VariantPredicates.ref_length('>', 5)
        >>> predicate.get_question()
        'ref allele length > 5'

        Args:
            operator: a `str` with the desired test. Must be one of ``{ '<', '<=', '==', '!=', '>=', '>' }``.
            length: a non-negative `int` with the length threshold.
        """
        return RefAlleleLengthPredicate(operator, length)

    @staticmethod
    def change_length(
        operator: typing.Literal["<", "<=", "==", "!=", ">=", ">"],
        threshold: int,
    ) -> VariantPredicate:
        """
        Prepare a :class:`VariantPredicate` for testing if the variant's change length
        is above, below, or (not) equal to certain `threshold`.

        .. seealso::

            See :ref:`change-length-of-an-allele` for more info.

        **Example**

        Make a predicate for testing if the change length is less than or equal to `-10`,
        e.g. to test if a variant is a *deletion* leading to removal of at least 10 base pairs:

        >>> from gpsea.analysis.predicate.genotype import VariantPredicates
        >>> predicate = VariantPredicates.change_length('<=', -10)
        >>> predicate.get_question()
        'change length <= -10'

        Args:
            operator: a `str` with the desired test. Must be one of ``{ '<', '<=', '==', '!=', '>=', '>' }``.
            threshold: an `int` with the threshold. Can be negative, zero, or positive.
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

            The predicate uses :meth:`~gpsea.model.VariantCoordinates.change_length`
            to determine if the length of the variant is above or below `threshold`.

            **IMPORTANT**: the change lengths of deletions are *negative*, since the alternate allele
            is shorter than the reference allele. See :ref:`change-length-of-an-allele` for more info.

        **Example**

        Prepare a predicate for testing if the variant is a chromosomal deletion that removes at least 20 base pairs:

        >>> from gpsea.analysis.predicate.genotype import VariantPredicates
        >>> predicate = VariantPredicates.is_structural_deletion(-20)
        >>> predicate.get_question()
        '(structural type is SO:1000029 OR (variant class is DEL AND change length <= -20))'

        Args:
            threshold: an `int` with the change length threshold to determine if a variant is "structural"
              (-50 bp by default).
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
        self, feature_type: FeatureType, tx_id: str
    ) -> VariantPredicate:
        """
        Prepare a :class:`VariantPredicate` that tests if the variant affects a protein feature type.

        Args:
            feature_type: the target protein :class:`FeatureType` (e.g. :class:`FeatureType.DOMAIN`)
        """
        return ProteinFeatureTypePredicate(
            feature_type, tx_id, self._protein_metadata_service
        )

    def protein_feature(self, feature_id: str, tx_id: str) -> VariantPredicate:
        """
        Prepare a :class:`VariantPredicate` that tests if the variant affects a protein feature type.

        Args:
            feature_id: the id of the target protein feature (e.g. `ANK 1`)
            tx_id: a `str` with the accession ID of the target transcript (e.g. `NM_123.4`)
        """
        return ProteinFeaturePredicate(
            feature_id, tx_id, self._protein_metadata_service
        )
