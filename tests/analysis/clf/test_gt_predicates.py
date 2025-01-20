import pytest

from gpsea.model import Patient, Sex, SampleLabels, VariantEffect
from gpsea.analysis.clf import (
    sex_classifier,
    monoallelic_classifier,
    biallelic_classifier,
    allele_count,
)
from gpsea.analysis.predicate import variant_effect


TX_ID = "tx:xyz"


class TestAlleleCount:
    @pytest.mark.parametrize(
        "patient_name,name",
        [
            ("adam", "1 allele"),
            ("eve", "1 allele"),
            ("cain", "1 allele"),
        ],
    )
    def test_ad_family__all_variants(
        self,
        patient_name: str,
        name: str,
        request: pytest.FixtureRequest,
    ):
        patient = request.getfixturevalue(patient_name)
        predicate = allele_count(counts=((0,), (1,)))

        categorization = predicate.test(patient)

        assert categorization is not None

        assert categorization.category.name == name

    @pytest.mark.parametrize(
        "patient_name,name",
        [
            ("adam", "0 alleles"),
            ("eve", "1 allele"),
            ("cain", "1 allele"),
        ],
    )
    def test_ad_family__missense_subset(
        self,
        patient_name: str,
        name: str,
        request: pytest.FixtureRequest,
    ):
        is_missense = variant_effect(VariantEffect.MISSENSE_VARIANT, tx_id=TX_ID)
        patient = request.getfixturevalue(patient_name)
        predicate = allele_count(
            counts=((0,), (1,)),
            target=is_missense,
        )

        categorization = predicate.test(patient)

        assert categorization is not None

        assert categorization.category.name == name

    @pytest.mark.parametrize(
        "patient_name,name",
        [
            ("walt", "2 alleles"),
            ("skyler", "2 alleles"),
            ("flynn", "2 alleles"),
            ("holly", "2 alleles"),
        ],
    )
    def test_ar_family__all_variants(
        self,
        patient_name: str,
        name: str,
        request: pytest.FixtureRequest,
    ):
        patient = request.getfixturevalue(patient_name)
        predicate = allele_count(counts=((0, 1), (2,)))

        categorization = predicate.test(patient)

        assert categorization is not None

        assert categorization.category.name == name

    @pytest.mark.parametrize(
        "patient_name,name",
        [
            ("walt", "1 allele"),
            ("skyler", "1 allele"),
            ("flynn", "2 alleles"),
            ("holly", "0 alleles"),
        ],
    )
    def test_ar_family__only_missense(
        self,
        patient_name: str,
        name: str,
        request: pytest.FixtureRequest,
    ):
        patient = request.getfixturevalue(patient_name)
        is_missense = variant_effect(VariantEffect.MISSENSE_VARIANT, tx_id=TX_ID)
        predicate = allele_count(
            counts=((0,), (1,), (2,)),
            target=is_missense,
        )

        categorization = predicate.test(patient)

        assert categorization is not None

        assert categorization.category.name == name

    def test_eq_and_hash(self):
        a = allele_count(counts=((0, 1), (2,)))
        b = allele_count(counts=((0, 1), (2,)))

        assert a == b
        assert hash(a) == hash(b)

    def test_summarize_groups(self):
        a = allele_count(counts=((0, 1), (2,)))

        assert a.summarize_classes() == "Allele count: 0 alleles OR 1 allele, 2 alleles"


class TestAllelePredicates:
    @pytest.mark.parametrize(
        "individual_name,expected_name",
        [
            ("adam", "B"),  # 0/0 & 0/1
            ("eve", "A"),  # 0/1 & 0/0
            ("cain", "A"),  # 0/1 & 0/0
        ],
    )
    def test_monoallelic_predicate_ad_family(
        self,
        individual_name: str,
        expected_name: str,
        request: pytest.FixtureRequest,
    ):
        is_missense = variant_effect(VariantEffect.MISSENSE_VARIANT, TX_ID)
        is_synonymous = variant_effect(VariantEffect.SYNONYMOUS_VARIANT, TX_ID)
        gt_clf = monoallelic_classifier(is_missense, is_synonymous)
        individual = request.getfixturevalue(individual_name)

        actual_cat = gt_clf.test(individual)

        assert actual_cat is not None
        assert actual_cat.category.name == expected_name

    def test_monoallelic_predicate__general_stuff(
        self,
    ):
        is_missense = variant_effect(VariantEffect.MISSENSE_VARIANT, TX_ID)
        is_synonymous = variant_effect(VariantEffect.SYNONYMOUS_VARIANT, TX_ID)

        gt_predicate = monoallelic_classifier(is_missense, is_synonymous)

        assert gt_predicate.summarize_classes() == "Allele group: A, B"

    @pytest.mark.parametrize(
        "individual_name,expected_name",
        [
            ("walt", "A/B"),  # 0/1 & 0/1
            ("skyler", "A/B"),  # 0/1 & 0/1
            ("flynn", "A/A"),  # 1/1 & 0/0
            ("holly", "B/B"),  # 0/0 & 1/1
        ],
    )
    def test_biallelic_predicate(
        self,
        individual_name: str,
        expected_name: str,
        request: pytest.FixtureRequest,
    ):
        is_missense = variant_effect(VariantEffect.MISSENSE_VARIANT, TX_ID)
        is_synonymous = variant_effect(VariantEffect.SYNONYMOUS_VARIANT, TX_ID)
        gt_predicate = biallelic_classifier(is_missense, is_synonymous)
        individual = request.getfixturevalue(individual_name)

        actual_cat = gt_predicate.test(individual)

        assert actual_cat is not None
        assert actual_cat.category.name == expected_name

    def test_biallelic_predicate__general_stuff(
        self,
    ):
        is_missense = variant_effect(VariantEffect.MISSENSE_VARIANT, TX_ID)
        is_synonymous = variant_effect(VariantEffect.SYNONYMOUS_VARIANT, TX_ID)

        gt_predicate = biallelic_classifier(is_missense, is_synonymous)

        assert gt_predicate.summarize_classes() == "Allele group: A/A, A/B, B/B"


class TestSexPredicate:
    def test_sex_predicate(
        self,
    ):
        joe = TestSexPredicate.make_patient("Joe", Sex.MALE)
        jane = TestSexPredicate.make_patient("Jane", Sex.FEMALE)
        miffy = TestSexPredicate.make_patient("Miffy", Sex.UNKNOWN_SEX)

        gt_predicate = sex_classifier()
        female, male = gt_predicate.get_categorizations()

        assert gt_predicate.test(joe) == male
        assert gt_predicate.test(jane) == female
        assert gt_predicate.test(miffy) is None

    def test_summarize_groups(self):
        gt_predicate = sex_classifier()

        assert gt_predicate.summarize_classes() == "Sex: FEMALE, MALE"

    @staticmethod
    def make_patient(label: str, sex: Sex) -> Patient:
        return Patient.from_raw_parts(
            SampleLabels(label),
            sex,
        )
