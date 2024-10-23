import pytest

import pandas as pd

from gpsea.analysis.predicate.genotype import GenotypePolyPredicate
from gpsea.analysis.pscore import PhenotypeScoreAnalysisResult
from gpsea.analysis.pscore.stats import MannWhitneyStatistic


class TestPhenotypeScoreAnalysisResult:

    @pytest.fixture(scope="class")
    def result(
        self,
        suox_gt_predicate: GenotypePolyPredicate,
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
            statistic=MannWhitneyStatistic(),
            data=data,
            pval=0.1234,
        )

    def test_properties(
        self,
        result: PhenotypeScoreAnalysisResult,
    ):
        assert tuple(result.gt_predicate.get_category_names()) == (
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
