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
            "Skip term if underpowered for 2x2 analysis",
            None,
            None,
        ]

    def test_mtc_filter_annotation_frequency_threshold_raises(
        self,
        hpo: hpotk.MinimalOntology,
    ):
        with pytest.raises(AssertionError) as e:
            IfHpoFilter.default_filter(
                hpo=hpo,
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
