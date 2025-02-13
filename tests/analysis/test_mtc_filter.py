import typing

import hpotk
import numpy as np
import pandas as pd
import pytest

from gpsea.analysis.mtc_filter import IfHpoFilter, SpecifiedTermsMtcFilter
from gpsea.analysis.clf import GenotypeClassifier, PhenotypeClassifier, HpoClassifier
from gpsea.analysis.pcats import apply_classifiers_on_individuals
from gpsea.model import Cohort


@pytest.fixture(scope="module")
def patient_counts(
    suox_cohort: Cohort,
    suox_gt_clf: GenotypeClassifier,
    suox_pheno_clfs: typing.Sequence[PhenotypeClassifier[hpotk.TermId]],
) -> typing.Sequence[pd.DataFrame]:
    _, counts = apply_classifiers_on_individuals(
        individuals=suox_cohort.all_patients,
        gt_clf=suox_gt_clf,
        pheno_clfs=suox_pheno_clfs,
    )
    return counts


class TestIfHpoFilter:
    @pytest.fixture
    def mtc_filter(
        self,
        hpo: hpotk.MinimalOntology,
    ) -> IfHpoFilter:
        return IfHpoFilter.default_filter(
            hpo=hpo,
            term_frequency_threshold=0.2,
            annotation_frequency_threshold=0.1,
        )

    @pytest.fixture(scope="class")
    def gt_clf(
        self,
        suox_gt_clf: GenotypeClassifier,
    ) -> GenotypeClassifier:
        return suox_gt_clf

    @pytest.fixture(scope="class")
    def ph_predicate(
        self,
        hpo: hpotk.MinimalOntology,
    ) -> PhenotypeClassifier[hpotk.TermId]:
        """
        For the purpose of testing counts, let's pretend the counts
        were created by this predicate.
        """
        return HpoClassifier(
            hpo=hpo,
            query=hpotk.TermId.from_curie("HP:0001250"),  # Seizure
            missing_implies_phenotype_excluded=False,
        )

    @staticmethod
    def prepare_counts_df(
        counts,
        gt_predicate: GenotypeClassifier,
        ph_predicate: PhenotypeClassifier[hpotk.TermId],
    ):
        values = np.array(counts).reshape((2, 2))
        index = pd.Index(ph_predicate.get_categories())
        columns = pd.Index(gt_predicate.get_categories())
        return pd.DataFrame(data=values, index=index, columns=columns)

    @pytest.mark.parametrize(
        "counts, expected",
        [
            ((10, 1, 15, 0), False),
            ((10, 0, 15, 1), False),
            ((10, 0, 15, 0), True),
            ((0, 15, 0, 19), True),
        ],
    )
    def test_one_genotype_has_zero_hpo_observations(
        self,
        counts: typing.Sequence[int],
        expected: bool,
        gt_clf: GenotypeClassifier,
        ph_predicate: PhenotypeClassifier[hpotk.TermId],
    ):
        counts_df = TestIfHpoFilter.prepare_counts_df(counts, gt_clf, ph_predicate)

        actual = IfHpoFilter.one_genotype_has_zero_hpo_observations(
            counts=counts_df,
            gt_clf=gt_clf,
        )

        assert actual == expected

    @pytest.mark.parametrize(
        "counts, expected",
        [
            ((1, 2, 99, 198), 0.01),
            ((1, 3, 99, 197), 0.015),
            ((0, 0, 100, 200), 0.0),
            ((0, 0, 0, 200), 0.0),
            ((0, 0, 0, 0), 0.0),
        ],
    )
    def test_get_maximum_group_observed_HPO_frequency(
        self,
        counts: typing.Tuple[int],
        expected: float,
        gt_clf: GenotypeClassifier,
        ph_predicate: PhenotypeClassifier[hpotk.TermId],
    ):
        counts_df = TestIfHpoFilter.prepare_counts_df(counts, gt_clf, ph_predicate)

        actual = IfHpoFilter.get_maximum_group_observed_HPO_frequency(
            counts_frame=counts_df,
            ph_clf=ph_predicate,
        )

        assert actual == pytest.approx(expected)

    def test_filter_terms_to_test(
        self,
        mtc_filter: IfHpoFilter,
        suox_gt_clf: GenotypeClassifier,
        suox_pheno_clfs: typing.Sequence[PhenotypeClassifier[hpotk.TermId]],
        patient_counts: typing.Sequence[pd.DataFrame],
        suox_cohort: Cohort,
    ):
        mtc_report = mtc_filter.filter(
            gt_clf=suox_gt_clf,
            pheno_clfs=suox_pheno_clfs,
            counts=patient_counts,
            cohort_size=len(suox_cohort),
        )

        assert isinstance(mtc_report, typing.Sequence)
        assert len(mtc_report) == 5

        is_ok = [r.is_passed() for r in mtc_report]
        assert is_ok == [True, True, False, True, True]

        reasons = [r.reason for r in mtc_report]
        assert reasons == [
            None,
            None,
            "Skipping term with less than 7 observations (not powered for 2x2)",
            None,
            None,
        ]

    def test_min_observed_HPO_threshold(
        self,
        suox_pheno_clfs: typing.Sequence[PhenotypeClassifier[hpotk.TermId]],
        patient_counts: typing.Sequence[pd.DataFrame],
    ):
        """
        In our heuristic filter, we only test terms that have at least a threshold
        frequency in at least one of the groups. We use the `patient_counts` - a sequence of DataFrames
        with 2x2 contingenicy tables of counts. For instance, each column will have one row for
        PatientCategories.YES and one for PatientCategories.NO, indicating counts of measured observed/excluded
        HPO phenotypes. Each column is a certain genotype, e.g., MISSENSE or NON-MISSENSE. We want the
        function to return the maximum frequency. In each column, the frequency is calculate by
        PatientCategories.YES / (PatientCategories.YES+PatientCategories.NO). This function tests that this works
        for all of the HPO terms in the dictionary.
        """
        EPSILON = 0.001
        curie2idx = {p.phenotype.value: i for i, p in enumerate(suox_pheno_clfs)}
        # Ectopia lentis HP:0001083  (1 2  3 1), freqs are 1/4=0.25 and 3/4=0.75
        idx = curie2idx["HP:0001083"]
        ectopia = patient_counts[idx]
        ectopia_predicate = suox_pheno_clfs[idx]
        max_f = IfHpoFilter.get_maximum_group_observed_HPO_frequency(
            ectopia,
            ph_clf=ectopia_predicate,
        )
        assert max_f == pytest.approx(0.75, abs=EPSILON)

        # Seizure HP:0001250 (11 5 0 1), freqs are 11/11=1.0 and 5/6=0.8333333
        idx = curie2idx["HP:0001250"]
        seizure = patient_counts[idx]
        seizure_predicate = suox_pheno_clfs[idx]
        max_f = IfHpoFilter.get_maximum_group_observed_HPO_frequency(
            seizure, ph_clf=seizure_predicate
        )
        assert max_f == pytest.approx(1.0, abs=EPSILON)

        # Sulfocysteinuria HP:0032350 (2 3 0 0), freqs are both 1
        idx = curie2idx["HP:0032350"]
        sulfocysteinuria = patient_counts[idx]
        sulfocysteinuria_predicate = suox_pheno_clfs[idx]
        max_f = IfHpoFilter.get_maximum_group_observed_HPO_frequency(
            sulfocysteinuria,
            ph_clf=sulfocysteinuria_predicate,
        )
        assert max_f == pytest.approx(1.0, abs=EPSILON)

        # Neurodevelopmental delay HP:0012758 (4 0 4 5), freqs are 4/8 = 0.5 and 0/5=0.0
        idx = curie2idx["HP:0012758"]
        ndelay = patient_counts[idx]
        ndelay_predicate = suox_pheno_clfs[idx]
        max_f = IfHpoFilter.get_maximum_group_observed_HPO_frequency(
            ndelay,
            ph_clf=ndelay_predicate,
        )
        assert max_f == pytest.approx(0.5, abs=EPSILON)

        # Hypertonia HP:0001276 (4 2 3 3) freqs are 4/7=0.4375 and 2/5=0.5714
        idx = curie2idx["HP:0001276"]
        hypertonia = patient_counts[idx]
        hypertonia_predicate = suox_pheno_clfs[idx]
        max_f = IfHpoFilter.get_maximum_group_observed_HPO_frequency(
            hypertonia,
            ph_clf=hypertonia_predicate,
        )
        assert max_f == pytest.approx(0.5714, abs=EPSILON)

    def test_mtc_filter_term_frequency_threshold_raises(
        self,
        hpo: hpotk.MinimalOntology,
    ):
        with pytest.raises(AssertionError) as e:
            IfHpoFilter.default_filter(
                hpo=hpo,
                term_frequency_threshold=1.1,
                annotation_frequency_threshold=0.1,
            )
        assert e.value.args == (
            "The term_frequency_threshold must be in the range (0, 1]",
        )

    def test_mtc_filter_annotation_frequency_threshold_raises(
        self,
        hpo: hpotk.MinimalOntology,
    ):
        with pytest.raises(AssertionError) as e:
            IfHpoFilter.default_filter(
                hpo=hpo,
                term_frequency_threshold=0.1,
                annotation_frequency_threshold=1.1,
            )
        assert e.value.args == (
            "The annotation_frequency_threshold must be in the range (0, 1]",
        )


class TestSpecifyTermsMtcFilter:
    def test_specified_term_mtc_filter(
        self,
        suox_gt_clf: GenotypeClassifier,
        suox_pheno_clfs: typing.Sequence[PhenotypeClassifier[hpotk.TermId]],
        patient_counts: typing.Sequence[pd.DataFrame],
        suox_cohort: Cohort,
    ):
        """
        The point of this test is to check that if we filter to test only one term ("HP:0032350"), then this
        is the only term that should survive the filter. We start with a total of five terms (n_usable==5),
        but after our filter, only one survives, and we have four cases in which the
        reason for filtering out is 'Non-specified term'
        """
        specified_filter = SpecifiedTermsMtcFilter(
            terms_to_test=(hpotk.TermId.from_curie("HP:0032350"),),
        )

        mtc_report = specified_filter.filter(
            gt_clf=suox_gt_clf,
            pheno_clfs=suox_pheno_clfs,
            counts=patient_counts,
            cohort_size=len(suox_cohort),
        )
        assert isinstance(mtc_report, typing.Sequence)
        assert len(mtc_report) == 5

        is_passed = [r.is_passed() for r in mtc_report]
        assert is_passed == [
            False,
            False,
            True,
            False,
            False,
        ]

        reasons = [r.reason for r in mtc_report]
        assert reasons == [
            "Non-specified term",
            "Non-specified term",
            None,
            "Non-specified term",
            "Non-specified term",
        ]

    @pytest.mark.parametrize(
        "val, msg",
        [
            ("NotACurie", "The CURIE NotACurie has no colon `:` or underscore `_`"),
            (0, "0 is neither `str` nor `hpotk.TermId`"),
        ],
    )
    def test_explodes_if_invalid_terms_provided(self, val: typing.Any, msg: str):
        with pytest.raises(ValueError) as e:
            SpecifiedTermsMtcFilter(
                terms_to_test=(val,),
            )
        assert e.value.args == (msg,)
