import typing

import pytest

from gpsea.model import *
from gpsea.analysis.predicate.genotype import (
    GenotypePolyPredicate,
    groups_predicate,
    filtering_predicate,
    VariantPredicates,
    VariantPredicate,
    ModeOfInheritancePredicate,
)


TX_ID = "tx:xyz"


class TestGroupsPredicate:

    @pytest.fixture(scope="class")
    def predicate(self) -> GenotypePolyPredicate:
        return groups_predicate(
            predicates=(
                VariantPredicates.variant_effect(VariantEffect.MISSENSE_VARIANT, TX_ID),
                VariantPredicates.variant_effect(
                    VariantEffect.FRAMESHIFT_VARIANT, TX_ID
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
        question = predicate.get_question_base()
        assert question == "Genotype group"

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
        patient_w_missense: Patient,
        predicate: GenotypePolyPredicate,
    ):
        cat = predicate.test(patient_w_missense)

        assert cat is not None
        assert cat.category.cat_id == 0
        assert cat.category.name == "Point"
        assert cat.category.description == "MISSENSE_VARIANT on tx:xyz"

    def test_test__frameshift(
        self,
        patient_w_frameshift: Patient,
        predicate: GenotypePolyPredicate,
    ):
        cat = predicate.test(patient_w_frameshift)

        assert cat is not None
        assert cat.category.cat_id == 1
        assert cat.category.name == "LoF"
        assert cat.category.description == "FRAMESHIFT_VARIANT on tx:xyz"


class TestModeOfInheritancePredicate:

    @pytest.fixture(scope="class")
    def variant_predicate(self) -> VariantPredicate:
        return VariantPredicates.variant_effect(VariantEffect.MISSENSE_VARIANT, TX_ID)

    @pytest.mark.parametrize(
        "patient_name,name",
        [
            ("adam", "HOM_REF"),
            ("eve", "HET"),
            ("cain", "HET"),
        ],
    )
    def test_autosomal_dominant(
        self,
        patient_name: str,
        name: str,
        variant_predicate: VariantPredicate,
        request: pytest.FixtureRequest,
    ):
        patient = request.getfixturevalue(patient_name)
        predicate = ModeOfInheritancePredicate.autosomal_dominant(variant_predicate)

        categorization = predicate.test(patient)

        assert categorization is not None

        assert categorization.category.name == name

    @pytest.mark.parametrize(
        "patient_name,name",
        [
            ("walt", "HET"),
            ("skyler", "HET"),
            ("flynn", "BIALLELIC_ALT"),
            ("holly", "HOM_REF"),
        ],
    )
    def test_autosomal_recessive(
        self,
        patient_name: str,
        name: str,
        variant_predicate: VariantPredicate,
        request: pytest.FixtureRequest,
    ):
        patient = request.getfixturevalue(patient_name)
        predicate = ModeOfInheritancePredicate.autosomal_recessive(variant_predicate)

        categorization = predicate.test(patient)

        assert categorization is not None

        assert categorization.category.name == name

    @pytest.mark.parametrize(
        "patient_name,name",
        [
            ("adam", "HOM_REF"),
            ("eve", "HET"),
            ("cain", "HET"),
        ],
    )
    def test_x_dominant(
        self,
        patient_name: str,
        name: str,
        variant_predicate: VariantPredicate,
        request: pytest.FixtureRequest,
    ):
        patient = request.getfixturevalue(patient_name)
        predicate = ModeOfInheritancePredicate.x_dominant(variant_predicate)

        categorization = predicate.test(patient)

        assert categorization is not None

        assert categorization.category.name == name

    @pytest.mark.parametrize(
        "patient_name,name",
        [
            ("anakin", "HOM_REF"),
            ("padme", "HET"),
            ("luke", "HEMI"),
            ("leia", "HET"),
        ],
    )
    def test_x_recessive(
        self,
        patient_name: str,
        name: str,
        variant_predicate: VariantPredicate,
        request: pytest.FixtureRequest,
    ):
        patient = request.getfixturevalue(patient_name)
        predicate = ModeOfInheritancePredicate.x_recessive(variant_predicate)

        categorization = predicate.test(patient)

        assert categorization is not None

        assert categorization.category.name == name


class TestPolyPredicate:

    @pytest.fixture(scope="class")
    def x_recessive_gt_predicate(self) -> GenotypePolyPredicate:
        affects_suox = VariantPredicates.gene("SUOX")
        return ModeOfInheritancePredicate.x_recessive(
            variant_predicate=affects_suox,
        )

    @pytest.mark.parametrize(
        "indices, expected",
        [
            ((0, 1), 2),
            ((1, 0), 2),
            ((1, 2), 2),
        ],
    )
    def test_filtering_predicate(
        self,
        indices: typing.Sequence[int],
        expected: int,
        x_recessive_gt_predicate: GenotypePolyPredicate,
    ):
        cats = x_recessive_gt_predicate.get_categorizations()
        targets = [cats[i] for i in indices]
        predicate = filtering_predicate(
            predicate=x_recessive_gt_predicate,
            targets=targets,
        )

        actual = len(predicate.get_categorizations())

        assert actual == expected

    def test_filtering_predicate__explodes_when_not_subsetting(
        self,
        x_recessive_gt_predicate: GenotypePolyPredicate,
    ):
        with pytest.raises(ValueError) as ve:
            filtering_predicate(
                predicate=x_recessive_gt_predicate,
                targets=x_recessive_gt_predicate.get_categorizations(),
            )

        assert (
            ve.value.args[0]
            == "It makes no sense to subset the a predicate with 4 categorizations with the same number (4) of targets"
        )

    def test_filtering_predicate__explodes_when_using_random_junk(
        self,
        x_recessive_gt_predicate: GenotypePolyPredicate,
    ):
        with pytest.raises(ValueError) as ve:
            filtering_predicate(
                predicate=x_recessive_gt_predicate,
                targets=(0, 1),
            )

        assert (
            ve.value.args[0]
            == "The targets at following indices are not categorizations: [0, 1]"
        )

    def test_filtering_predicate__explodes_when_using_one_category(
        self,
        x_recessive_gt_predicate: GenotypePolyPredicate,
    ):
        with pytest.raises(ValueError) as ve:
            filtering_predicate(
                predicate=x_recessive_gt_predicate,
                targets=(x_recessive_gt_predicate.get_categorizations()[0],),
            )

        assert (
            ve.value.args[0]
            == "At least 2 target categorizations must be provided but got 1"
        )
