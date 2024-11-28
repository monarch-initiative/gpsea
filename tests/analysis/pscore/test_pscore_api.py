import random
import pytest

import pandas as pd

from gpsea.analysis import StatisticResult
from gpsea.analysis.predicate.genotype import GenotypePolyPredicate
from gpsea.analysis.pscore import PhenotypeScoreAnalysisResult, PhenotypeScorer
from gpsea.analysis.pscore.stats import MannWhitneyStatistic


class TestPhenotypeScoreAnalysisResult:

    @pytest.fixture(scope="class")
    def phenotype_scorer(self) -> PhenotypeScorer:
        return PhenotypeScorer.wrap_scoring_function(
            func=lambda patient: random.random(),
            name="Random phenotype scorer",
        )

    @pytest.fixture(scope="class")
    def result(
        self,
        suox_gt_predicate: GenotypePolyPredicate,
        phenotype_scorer: PhenotypeScorer,
    ) -> PhenotypeScoreAnalysisResult:
        data = pd.DataFrame(
            data={
                "patient_id": ["A", "B", "C"],
                "genotype": [0, 1, None],
                "phenotype": [
                    10.,
                    float('nan'),
                    -4.,
                ],
            }
        ).set_index("patient_id")
        return PhenotypeScoreAnalysisResult(
            gt_predicate=suox_gt_predicate,
            phenotype=phenotype_scorer,
            statistic=MannWhitneyStatistic(),
            data=data,
            statistic_result=StatisticResult(statistic=.2, pval=0.1234),
        )

    def test_properties(
        self,
        result: PhenotypeScoreAnalysisResult,
    ):
        assert tuple(result.gt_predicate.group_labels) == (
            "0",
            "1",
        )
        assert result.data.index.to_list() == ["A", "B", "C"]
        assert result.pval == pytest.approx(0.1234)

    def test_complete_records(
        self,
        result: PhenotypeScoreAnalysisResult,
    ):
        records = result.complete_records()

        assert records.shape == (1, 2)
        assert records.loc["A", "genotype"] == 0
        assert records.loc["A", "phenotype"] == pytest.approx(10.)
