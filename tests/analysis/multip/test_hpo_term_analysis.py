import typing

import hpotk
import pytest

from gpsea.model import Cohort

from gpsea.analysis.mtc_filter import PhenotypeMtcFilter, HpoMtcFilter
from gpsea.analysis.multip import HpoTermAnalysis
from gpsea.analysis.predicate.genotype import GenotypePolyPredicate
from gpsea.analysis.predicate.phenotype import PhenotypePolyPredicate
from gpsea.analysis.stats import CountStatistic, ScipyFisherExact


class TestHpoTermAnalysis:

    @pytest.fixture(scope='class')
    def count_statistic(self) -> CountStatistic:
        return ScipyFisherExact()

    @pytest.fixture(scope='class')
    def phenotype_mtc_filter(
        self,
        hpo: hpotk.MinimalOntology,
    ) -> PhenotypeMtcFilter:
        return HpoMtcFilter.default_filter(
            hpo=hpo,
            term_frequency_threshold=.2,
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
            mtc_correction='fdr_bh',
            mtc_alpha=.05,
        )

    def test_compare_genotype_vs_phenotypes(
        self,
        analysis: HpoTermAnalysis,
        suox_cohort: Cohort,
        suox_gt_predicate: GenotypePolyPredicate,
        suox_pheno_predicates: typing.Sequence[PhenotypePolyPredicate[hpotk.TermId]]
    ):
        results = analysis.compare_genotype_vs_phenotypes(
            cohort=suox_cohort.all_patients,
            gt_predicate=suox_gt_predicate,
            pheno_predicates=suox_pheno_predicates,
        )
        
        assert results is not None
        # TODO: improve testing
