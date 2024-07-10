import pytest

from genophenocorr.analysis.predicate.genotype import VariantPredicates, ProteinPredicates
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


class TestBuildGenotypePredicate:

    def test_build_it(self):
        """
        """
        predicate = VariantPredicates.gene(symbol='A').oder(VariantPredicates.gene(symbol='B'))
        print(predicate)

