import pytest

from gpsea.model import *
from gpsea.analysis.predicate.genotype import (
    GenotypePolyPredicate,
    groups_predicate,
    sex_predicate,
    monoallelic_predicate,
    biallelic_predicate,
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


class TestAllelePredicates:

    @pytest.mark.parametrize(
        "individual_name,expected_name",
        [
            ("adam", "Zero"),  # 0/0
            ("eve", "One"),  # 0/1
            ("cain", "One"),  # 0/1
        ],
    )
    def test_monoallelic_predicate_ad_family(
        self,
        individual_name: str,
        expected_name: str,
        request: pytest.FixtureRequest,
    ):
        is_missense = VariantPredicates.variant_effect(VariantEffect.MISSENSE_VARIANT, TX_ID)
        gt_predicate = monoallelic_predicate(is_missense)
        individual = request.getfixturevalue(individual_name)

        actual_cat = gt_predicate.test(individual)

        assert actual_cat is not None
        assert actual_cat.category.name == expected_name

    @pytest.mark.parametrize(
        "individual_name,not_none,expected_name",
        [
            ("walt", True, "One"),  # 1/2
            ("skyler", True, "One"),  # 1/2
            ("flynn", False, None),  # 1/1, hence 2 alleles which is too much to work with monoallelic predicate.
            ("holly", True, "Zero"),  # 2/2
        ],
    )
    def test_monoallelic_predicate_ar_family(
        self,
        individual_name: str,
        not_none: bool,
        expected_name: str,
        request: pytest.FixtureRequest,
    ):
        is_missense = VariantPredicates.variant_effect(VariantEffect.MISSENSE_VARIANT, TX_ID)
        gt_predicate = monoallelic_predicate(is_missense)
        individual = request.getfixturevalue(individual_name)

        actual_cat = gt_predicate.test(individual)

        if not_none:
            assert actual_cat is not None
            assert actual_cat.category.name == expected_name
        else:
            assert actual_cat is None

    @pytest.mark.parametrize(
        "individual_name,expected_name",
        [
            ("walt", "AB"),  # 1/2
            ("skyler", "AB"),  # 1/2
            ("flynn", "AA"),  # 1/1
            ("holly", "BB"),  # 2/2
        ],
    )
    def test_biallelic_predicate(
        self,
        individual_name: str,
        expected_name: str,
        request: pytest.FixtureRequest,
    ):
        is_missense = VariantPredicates.variant_effect(VariantEffect.MISSENSE_VARIANT, TX_ID)
        is_stop_gained = VariantPredicates.variant_effect(VariantEffect.STOP_GAINED, TX_ID)
        gt_predicate = biallelic_predicate(is_missense, is_stop_gained)
        individual = request.getfixturevalue(individual_name)

        actual_cat = gt_predicate.test(individual)

        assert actual_cat is not None
        assert actual_cat.category.name == expected_name


class TestSexPredicate:

    def test_sex_predicate(
        self,
    ):
        joe = TestSexPredicate.make_patient('Joe', Sex.MALE)
        jane = TestSexPredicate.make_patient('Jane', Sex.FEMALE)
        miffy = TestSexPredicate.make_patient('Miffy', Sex.UNKNOWN_SEX)
        
        gt_predicate = sex_predicate()
        female, male = gt_predicate.get_categorizations()

        assert gt_predicate.test(joe) == male
        assert gt_predicate.test(jane) == female
        assert gt_predicate.test(miffy) is None

    def test_get_question(self):
        gt_predicate = sex_predicate()
        assert gt_predicate.display_question() == 'Sex of the individual: FEMALE, MALE'

    @staticmethod
    def make_patient(label: str, sex: Sex) -> Patient:
        return Patient(
            SampleLabels(label),
            sex,
            phenotypes=(),
            diseases=(),
            variants=(),
        )
