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

    def test_is_large_imprecise_sv(
        self,
        variant: Variant,
        structural_variant: Variant,
    ):
        predicate = VariantPredicates.is_large_imprecise_sv()

        assert predicate.test(variant) == False
        assert predicate.test(structural_variant) == True

    def test_is_structural_predicate(
        self,
        variant: Variant,
        structural_variant: Variant,
    ):
        predicate = VariantPredicates.is_structural_variant()

        assert predicate.test(variant) == False
        assert predicate.test(structural_variant) == True

    def test_structural_type(
        self,
        variant: Variant,
        structural_variant: Variant,
    ):
        predicate = VariantPredicates.structural_type('SO:1000029')

        assert predicate.test(variant) == False
        assert predicate.test(structural_variant) == True

    def test_variant_class(
        self,
        variant: Variant,
    ):
        predicate = VariantPredicates.variant_class(VariantClass.SNV)

        assert predicate.test(variant) == True

    def test_change_length(
        self,
        variant: Variant,
        structural_variant: Variant,
    ):
        predicate = VariantPredicates.change_length('==', 0)
        
        # variant is an SNP
        assert predicate.test(variant) == True

        # structural_variant is an imprecise DEL, hence False
        assert predicate.test(structural_variant) == False

    def test_change_length_is_false_for_imprecise_SVs_no_matter_what(
        self,
        structural_variant: Variant,
    ):
        for op in ("<", "<=", "==", "!=", ">=", ">"):
            predicate = VariantPredicates.change_length(op, 0)
            assert predicate.test(structural_variant) == False

    def test_structural_deletion(
        self,
        structural_variant: Variant,
    ):
        predicate = VariantPredicates.is_structural_deletion()

        assert predicate.test(structural_variant) == True


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

    @pytest.mark.parametrize(
        'tx_id,expected',
        [
            ('tx:abc', False),
            ('whatever', True),
        ]
    )
    def test_inv_predicate(
        self,
        variant: Variant,
        tx_id: str,
        expected: bool,
    ):
        predicate = ~VariantPredicates.transcript(tx_id)

        assert predicate.test(variant) == expected
    
    def test_no_double_inv_happens(
        self,
    ):
        predicate = VariantPredicates.gene('FBN1')
        
        # Inverting a predicate must produce a new predicate.
        inv_predicate = ~predicate
        assert inv_predicate is not predicate

        # No double negation!
        inv_inv_predicate = ~inv_predicate
        assert inv_inv_predicate is predicate

    def test_empty_all_predicate_raises_error(
        self,
    ):
        with pytest.raises(ValueError) as e:
            empty = ()
            VariantPredicates.all(empty)
        assert e.value.args[0] == 'Predicates must not be empty!'
    
    def test_empty_any_predicate_raises_error(
        self,
    ):
        with pytest.raises(ValueError) as e:
            empty = ()
            VariantPredicates.any(empty)
        assert e.value.args[0] == 'Predicates must not be empty!'
