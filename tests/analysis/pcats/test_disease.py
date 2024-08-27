import hpotk
import pytest

from gpsea.model import Cohort

from gpsea.analysis.pcats import DiseaseAnalysis
from gpsea.analysis.pcats.stats import CountStatistic, ScipyFisherExact
from gpsea.analysis.predicate.genotype import GenotypePolyPredicate
from gpsea.analysis.predicate.phenotype import DiseasePresencePredicate


class TestDiseaseAnalysis:

    @pytest.fixture(scope='class')
    def count_statistic(self) -> CountStatistic:
        return ScipyFisherExact()

    @pytest.fixture(scope='class')
    def suox_disease(self) -> DiseasePresencePredicate:
        sulfite_oxidase_deficiency = hpotk.TermId.from_curie('OMIM:272300')
        return DiseasePresencePredicate(
            disease_id_query=sulfite_oxidase_deficiency,
        )

    @pytest.fixture
    def analysis(
        self,
        count_statistic: CountStatistic,
    ) -> DiseaseAnalysis:
        return DiseaseAnalysis(
            count_statistic=count_statistic,
            mtc_correction='fdr_bh',
            mtc_alpha=.05,
        )

    def test_compare_genotype_vs_phenotypes(
        self,
        analysis: DiseaseAnalysis,
        suox_cohort: Cohort,
        suox_gt_predicate: GenotypePolyPredicate,
        suox_disease: DiseasePresencePredicate,
    ):
        results = analysis.compare_genotype_vs_phenotypes(
            cohort=suox_cohort.all_patients,
            gt_predicate=suox_gt_predicate,
            pheno_predicates=(suox_disease,),
        )

        assert results is not None
        # TODO: improve testing
