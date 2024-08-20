import pytest

from genophenocorr.model import *
from genophenocorr.analysis.predicate import GenotypePolyPredicate
from genophenocorr.analysis.predicate.genotype import (
    groups_predicate,
    VariantPredicates,
)


class TestGroupsPredicate:

    TX_ID = "tx:xyz"

    @pytest.fixture(scope="class")
    def predicate(self) -> GenotypePolyPredicate:
        return groups_predicate(
            predicates=(
                VariantPredicates.variant_effect(
                    VariantEffect.MISSENSE_VARIANT, TestGroupsPredicate.TX_ID
                ),
                VariantPredicates.variant_effect(
                    VariantEffect.FRAMESHIFT_VARIANT, TestGroupsPredicate.TX_ID
                ),
            ),
            group_names=(
                "Point",
                "LoF",
            ),
        )

    def test_get_question(
        self,
        predicate: GenotypePolyPredicate,
    ):
        question = predicate.get_question()
        assert question == "What group does the patient belong to: Point, LoF"

    def test_get_categorizations(
        self,
        predicate: GenotypePolyPredicate,
    ):
        categorizations = predicate.get_categorizations()

        names = [c.category.name for c in categorizations]
        assert names == ["Point", "LoF"]

        descriptions = [c.category.description for c in categorizations]
        assert descriptions == [
            "MISSENSE_VARIANT on tx:xyz",
            "FRAMESHIFT_VARIANT on tx:xyz",
        ]

    def test_test__missense(
        self,
        patient_w_missense: Variant,
        predicate: GenotypePolyPredicate,
    ):
        cat = predicate.test(patient_w_missense)

        assert cat is not None
        assert cat.category.cat_id == 0
        assert cat.category.name == "Point"
        assert cat.category.description == "MISSENSE_VARIANT on tx:xyz"

    def test_test__frameshift(
        self,
        patient_w_frameshift: Variant,
        predicate: GenotypePolyPredicate,
    ):
        cat = predicate.test(patient_w_frameshift)

        assert cat is not None
        assert cat.category.cat_id == 1
        assert cat.category.name == "LoF"
        assert cat.category.description == "FRAMESHIFT_VARIANT on tx:xyz"
