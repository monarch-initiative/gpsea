import json
import os

import hpotk
import pytest

import pandas as pd

from gpsea.model import Cohort
from gpsea.io import GpseaJSONDecoder
from gpsea.analysis.predicate.genotype import (
    GenotypePolyPredicate,
    VariantPredicates,
    monoallelic_predicate,
)
from gpsea.analysis.temporal import SurvivalAnalysis, SurvivalAnalysisResult, Survival
from gpsea.analysis.temporal.endpoint import hpo_onset
from gpsea.analysis.temporal.stats import LogRankTest


@pytest.fixture(scope="module")
def umod_cohort(
    fpath_cohort_data_dir: str,
) -> Cohort:
    fpath_cohort = os.path.join(fpath_cohort_data_dir, "UMOD.0.1.20.json")
    with open(fpath_cohort) as fh:
        return json.load(fh, cls=GpseaJSONDecoder)


@pytest.fixture(scope="module")
def umod_gt_predicate() -> GenotypePolyPredicate:
    in_exon_3 = VariantPredicates.exon(3, tx_id="NM_003361.4")
    return monoallelic_predicate(
        a_predicate=in_exon_3, b_predicate=~in_exon_3,
        a_label="Exon 3", b_label="Other exon",
    )


class TestSurvivalAnalysis:

    @pytest.fixture(scope="class")
    def survival_analysis(self) -> SurvivalAnalysis:
        return SurvivalAnalysis(statistic=LogRankTest())

    def test_compare_genotype_vs_survival(
        self,
        survival_analysis: SurvivalAnalysis,
        hpo: hpotk.MinimalOntology,
        umod_cohort: Cohort,
        umod_gt_predicate: GenotypePolyPredicate,
    ):

        endpoint = hpo_onset(
            hpo=hpo,
            term_id="HP:0003774",  # Stage 5 chronic kidney disease
        )

        result = survival_analysis.compare_genotype_vs_survival(
            cohort=umod_cohort,
            gt_predicate=umod_gt_predicate,
            endpoint=endpoint,
        )

        assert result.pval == pytest.approx(0.062004258)


class TestSurvivalAnalysisResult:

    @pytest.fixture(scope="class")
    def result(
        self,
        umod_gt_predicate: GenotypePolyPredicate,
    ) -> SurvivalAnalysisResult:
        data = pd.DataFrame(
            data={
                "patient_id": ["A", "B", "C"],
                "genotype": [0, 1, None],
                "phenotype": [
                    Survival(value=12.3, is_censored=False),
                    None,
                    Survival(value=45.6, is_censored=True),
                ],
            }
        ).set_index("patient_id")
        return SurvivalAnalysisResult(
            gt_predicate=umod_gt_predicate,
            statistic=LogRankTest(),
            data=data,
            pval=0.1234,
        )

    def test_properties(
        self,
        result: SurvivalAnalysisResult,
    ):
        assert tuple(result.gt_predicate.get_category_names()) == (
            "Exon 3",
            "Other exon",
        )
        assert result.data.index.to_list() == ["A", "B", "C"]
        assert result.pval == pytest.approx(0.1234)

    def test_complete_records(
        self,
        result: SurvivalAnalysisResult,
    ):
        records = result.complete_records()

        assert records.shape == (1, 2)
        assert records.loc["A", "genotype"] == 0
        assert records.loc["A", "phenotype"] == Survival(value=12.3, is_censored=False)
