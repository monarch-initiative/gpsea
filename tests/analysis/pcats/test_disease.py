import hpotk
import pytest

from gpsea.model import Cohort

from gpsea.analysis.pcats import DiseaseAnalysis
from gpsea.analysis.pcats.stats import CountStatistic, FisherExactTest
from gpsea.analysis.clf import GenotypeClassifier, DiseasePresenceClassifier


class TestDiseaseAnalysis:
    @pytest.fixture(scope="class")
    def count_statistic(self) -> CountStatistic:
        return FisherExactTest()

    @pytest.fixture(scope="class")
    def suox_disease(self) -> DiseasePresenceClassifier:
        sulfite_oxidase_deficiency = hpotk.TermId.from_curie("OMIM:272300")
        return DiseasePresenceClassifier(
            disease_id_query=sulfite_oxidase_deficiency,
        )

    @pytest.fixture
    def analysis(
        self,
        count_statistic: CountStatistic,
    ) -> DiseaseAnalysis:
        return DiseaseAnalysis(
            count_statistic=count_statistic,
            mtc_correction="fdr_bh",
            mtc_alpha=0.05,
        )

    def test_compare_genotype_vs_phenotypes(
        self,
        analysis: DiseaseAnalysis,
        suox_cohort: Cohort,
        suox_gt_clf: GenotypeClassifier,
        suox_disease: DiseasePresenceClassifier,
    ):
        results = analysis.compare_genotype_vs_phenotypes(
            cohort=suox_cohort.all_patients,
            gt_clf=suox_gt_clf,
            pheno_clfs=(suox_disease,),
        )

        assert results is not None
        # TODO: improve testing
