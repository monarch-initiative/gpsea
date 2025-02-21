import typing

import hpotk
import numpy as np
import pytest

from gpsea.model import Cohort

from gpsea.analysis.mtc_filter import PhenotypeMtcFilter, IfHpoFilter
from gpsea.analysis.pcats import HpoTermAnalysis
from gpsea.analysis.pcats.stats import CountStatistic, FisherExactTest
from gpsea.analysis.clf import GenotypeClassifier, PhenotypeClassifier


class TestHpoTermAnalysis:
    @pytest.fixture(scope="class")
    def count_statistic(self) -> CountStatistic:
        return FisherExactTest()

    @pytest.fixture(scope="class")
    def phenotype_mtc_filter(
        self,
        hpo: hpotk.MinimalOntology,
    ) -> PhenotypeMtcFilter:
        return IfHpoFilter.default_filter(
            hpo=hpo,
            annotation_frequency_threshold=0.25,
        )

    @pytest.fixture
    def analysis(
        self,
        count_statistic: CountStatistic,
        phenotype_mtc_filter: PhenotypeMtcFilter,
    ) -> HpoTermAnalysis:
        return HpoTermAnalysis(
            count_statistic=count_statistic,
            mtc_filter=phenotype_mtc_filter,
            mtc_correction="fdr_bh",
            mtc_alpha=0.05,
        )

    def test_compare_genotype_vs_phenotypes(
        self,
        analysis: HpoTermAnalysis,
        suox_cohort: Cohort,
        suox_gt_clf: GenotypeClassifier,
        suox_pheno_clfs: typing.Sequence[PhenotypeClassifier[hpotk.TermId]],
    ):
        results = analysis.compare_genotype_vs_phenotypes(
            cohort=suox_cohort.all_patients,
            gt_clf=suox_gt_clf,
            pheno_clfs=suox_pheno_clfs,
        )

        assert results is not None

        assert results.total_tests == 3
        assert results.n_usable == (17, 7, 5, 13, 12)
        assert results.pvals == pytest.approx(
            [
                0.35294117647058815,
                float("nan"),
                float("nan"),
                0.1048951048951049,
                1.0,
            ],
            nan_ok=True,
        )
        assert results.corrected_pvals is not None
        assert results.corrected_pvals == pytest.approx(
            [
                0.5294117647058822,
                float("nan"),
                float("nan"),
                0.3146853146853147,
                1.0,
            ],
            nan_ok=True,
        )

    def test_compare_genotype_vs_phenotypes_can_handle_if_no_phenotypes_are_left_after_mtc_filter(
        self,
        analysis: HpoTermAnalysis,
        degenerated_cohort: Cohort,
        suox_gt_clf: GenotypeClassifier,
        suox_pheno_clfs: typing.Sequence[PhenotypeClassifier[hpotk.TermId]],
    ):
        result = analysis.compare_genotype_vs_phenotypes(
            cohort=degenerated_cohort,
            gt_clf=suox_gt_clf,
            pheno_clfs=suox_pheno_clfs,
        )

        assert (
            result.total_tests == 0
        ), "No tests should have been done due to MTC filtering"
        assert np.all(np.isnan(result.pvals)), "All p values should be NaN"
        assert result.corrected_pvals is None
