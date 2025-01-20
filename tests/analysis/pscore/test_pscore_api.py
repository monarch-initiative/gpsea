import random
import pytest

import pandas as pd

from gpsea.analysis import StatisticResult
from gpsea.analysis.clf import GenotypeClassifier
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
        suox_gt_clf: GenotypeClassifier,
        phenotype_scorer: PhenotypeScorer,
    ) -> PhenotypeScoreAnalysisResult:
        data = pd.DataFrame(
            data={
                "patient_id": [
                    "A",
                    "B",
                    "C",
                    "D",
                    "E",
                    "F",
                    "G",
                    "H",
                    "I",
                    "J",
                    "K",
                    "L",
                    "M",
                    "N",
                ],
                "genotype": [0, 1, None, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0],
                "phenotype": [
                    10.0,
                    float("nan"),
                    -4.0,
                    15.0,
                    float("nan"),
                    11.0,
                    2,
                    7,
                    -3,
                    16.0,
                    14.0,
                    9.0,
                    4.0,
                    6.0,
                ],
            }
        ).set_index("patient_id")
        return PhenotypeScoreAnalysisResult(
            gt_clf=suox_gt_clf,
            phenotype=phenotype_scorer,
            statistic=MannWhitneyStatistic(),
            data=data,
            statistic_result=StatisticResult(statistic=0.2, pval=0.1234),
        )

    def test_properties(
        self,
        result: PhenotypeScoreAnalysisResult,
    ):
        assert tuple(result.gt_clf.class_labels) == (
            "0 alleles",
            "1 allele",
        )
        assert result.data.index.to_list() == [
            "A",
            "B",
            "C",
            "D",
            "E",
            "F",
            "G",
            "H",
            "I",
            "J",
            "K",
            "L",
            "M",
            "N",
        ]
        assert result.pval == pytest.approx(0.1234)

    def test_complete_records(
        self,
        result: PhenotypeScoreAnalysisResult,
    ):
        records = result.complete_records()

        assert records.shape == (11, 2)
        assert records.loc["A", "genotype"] == 0
        assert records.loc["A", "phenotype"] == pytest.approx(10.0)

    @pytest.mark.skip("Run manually")
    def test_plot_boxplots(
        self,
        result: PhenotypeScoreAnalysisResult,
    ):
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(10, 6))

        result.plot_boxplots(
            ax=ax,
        )
        fig.savefig("boxplot.png")

    @pytest.mark.skip("Run manually")
    def test_plot_violins(
        self,
        result: PhenotypeScoreAnalysisResult,
    ):
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(10, 6))

        result.plot_violins(
            ax=ax,
        )
        fig.savefig("violinplot.png")
