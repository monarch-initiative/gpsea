import typing

import hpotk
import pytest

from gpsea.model import Cohort

from gpsea.analysis.mtc_filter import PhenotypeMtcFilter, HpoMtcFilter
from gpsea.analysis.pcats import HpoTermAnalysis
from gpsea.analysis.pcats.stats import CountStatistic, FisherExactTest
from gpsea.analysis.predicate.genotype import GenotypePolyPredicate
from gpsea.analysis.predicate.phenotype import PhenotypePolyPredicate


class TestHpoTermAnalysis:

    @pytest.fixture(scope="class")
    def count_statistic(self) -> CountStatistic:
        return FisherExactTest()

    @pytest.fixture(scope="class")
    def phenotype_mtc_filter(
        self,
        hpo: hpotk.MinimalOntology,
    ) -> PhenotypeMtcFilter:
        return HpoMtcFilter.default_filter(
            hpo=hpo,
            term_frequency_threshold=0.2,
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
        suox_gt_predicate: GenotypePolyPredicate,
        suox_pheno_predicates: typing.Sequence[PhenotypePolyPredicate[hpotk.TermId]],
    ):
        results = analysis.compare_genotype_vs_phenotypes(
            cohort=suox_cohort.all_patients,
            gt_predicate=suox_gt_predicate,
            pheno_predicates=suox_pheno_predicates,
        )

        assert results is not None

        assert results.total_tests == 4
        assert results.n_usable == (17, 7, 5, 13, 12)
        assert results.pvals == pytest.approx(
            [
                0.35294117647058815,
                0.48571428571428565,
                float("nan"),
                0.1048951048951049,
                1.,
            ],
            nan_ok=True,
        )
        assert results.corrected_pvals is not None
        assert results.corrected_pvals == pytest.approx(
            [
                0.6476190476190475,
                0.6476190476190475,
                float("nan"),
                0.4195804195804196,
                1.0,
            ],
            nan_ok=True,
        )

    def test_compare_genotype_vs_phenotypes_explodes_if_no_phenotypes_are_left_after_mtc_filter(
        self,
        analysis: HpoTermAnalysis,
        degenerated_cohort: Cohort,
        suox_gt_predicate: GenotypePolyPredicate,
        suox_pheno_predicates: typing.Sequence[PhenotypePolyPredicate[hpotk.TermId]],
    ):
        with pytest.raises(ValueError) as e:
            analysis.compare_genotype_vs_phenotypes(
                cohort=degenerated_cohort,
                gt_predicate=suox_gt_predicate,
                pheno_predicates=suox_pheno_predicates,
            )

        assert (
            e.value.args[0]
            == "No phenotypes are left for the analysis after MTC filtering step"
        )
