import pytest

from genophenocorr.analysis.predicate.genotype import VariantPredicates, ProteinPredicates, VariantPredicate
from genophenocorr.model import *
from genophenocorr.model.genome import *


class TestVariantPredicates:

    @pytest.mark.parametrize(
        'effect, expected',
        [
            (VariantEffect.MISSENSE_VARIANT, True),
            (VariantEffect.SPLICE_DONOR_VARIANT, True),
            (VariantEffect.STOP_GAINED, False),
            (VariantEffect.STOP_LOST, False),
        ]
    )
    def test_variant_effect_predicate(
            self,
            variant: Variant,
            effect: VariantEffect,
            expected: bool,
    ):
        predicate = VariantPredicates.variant_effect(effect, tx_id='tx:xyz')

        assert predicate.test(variant) == expected

    @pytest.mark.parametrize(
        'variant_key, expected',
        [
            ('22_101_101_C_G', True),
            ('22_101_101_C_T', False),
        ]
    )
    def test_variant_key_predicate(
            self,
            variant: Variant,
            variant_key: str,
            expected: bool,
    ):
        predicate = VariantPredicates.variant_key(variant_key)

        assert predicate.test(variant) == expected

    @pytest.mark.parametrize(
        'exon, expected',
        [
            (0, False),
            (4, True),
            (5, False),
        ]
    )
    def test_exon_predicate(
            self,
            variant: Variant,
            exon: int,
            expected: bool,
    ):
        predicate = VariantPredicates.exon(exon, tx_id='tx:xyz')

        assert predicate.test(variant) == expected

    @pytest.mark.parametrize(
        'tx_id, expected',
        [
            ('tx:xyz', True),
            ('other', False),
        ]
    )
    def test_transcript_predicate(
            self,
            variant: Variant,
            tx_id: str,
            expected: bool,
    ):
        predicate = VariantPredicates.transcript(tx_id)

        assert predicate.test(variant) == expected

    @pytest.mark.parametrize(
        'symbol, expected',
        [
            ('a_gene', True),
            ('b_gene', False),
        ]
    )
    def test_gene_predicate(
            self,
            variant: Variant,
            symbol: str,
            expected: bool,
    ):
        predicate = VariantPredicates.gene(symbol)

        assert predicate.test(variant) == expected


class TestProteinPredicates:

    @pytest.mark.parametrize(
        'feature_type, expected',
        [
            (FeatureType.DOMAIN, True),
            (FeatureType.REPEAT, False),
        ]
    )
    def test_protein_feature_type(
            self,
            protein_predicates: ProteinPredicates,
            variant: Variant,
            feature_type: FeatureType,
            expected: bool,
    ):
        predicate = protein_predicates.protein_feature_type(feature_type, tx_id='tx:xyz')

        assert predicate.test(variant) == expected

    @pytest.mark.parametrize(
        'feature_id, expected',
        [
            ('MOCK_DOMAIN', True),
            ('REAL_DOMAIN', False),
            ('MOCK_REPEAT', False),
        ]
    )
    def test_protein_feature_id(
            self,
            protein_predicates: ProteinPredicates,
            variant: Variant,
            feature_id: str,
            expected: bool,
    ):
        predicate = protein_predicates.protein_feature(feature_id, tx_id='tx:xyz')

        assert predicate.test(variant) == expected

class TestLogicalVariantPredicate:
    """
    Test that the AND and OR variant predicate combinators work as expected.
    """

    def test_equivalent_predicates_are_not_chained(self):
        a1 = VariantPredicates.gene(symbol='A')
        a2 = VariantPredicates.gene(symbol='A')

        assert a1 & a2 is a1
        assert a1 | a2 is a1

        assert a2 & a1 is a2
        assert a2 | a1 is a2

    @pytest.mark.parametrize(
        'left,right,expected',
        [
            ('tx:abc', 'tx:xyz', True),
            ('tx:abc', 'whatever', False),
            ('whatever', 'tx:xyz', False),
            ('whatever', 'whoever', False),
        ]
    )
    def test_und_predicate(
        self,
        variant: Variant,
        left: str,
        right: str,
        expected: bool,
    ):
        predicate = VariantPredicates.transcript(tx_id=left) & VariantPredicates.transcript(tx_id=right)

        assert predicate.test(variant) == expected

    @pytest.mark.parametrize(
        'left,right,expected',
        [
            ('tx:abc', 'tx:xyz', True),
            ('tx:abc', 'whatever', True),
            ('whatever', 'tx:xyz', True),
            ('whatever', 'whoever', False),
        ]
    )
    def test_or_predicate(
        self,
        variant: Variant,
        left: str,
        right: str,
        expected: bool,
    ):
        predicate = VariantPredicates.transcript(tx_id=left) | VariantPredicates.transcript(tx_id=right)
        
        assert predicate.test(variant) == expected
