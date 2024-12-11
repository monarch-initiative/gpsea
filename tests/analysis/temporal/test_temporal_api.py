import json
import os

import hpotk
import pytest

import matplotlib.pyplot as plt
import pandas as pd

from gpsea.analysis._base import AnalysisException
from gpsea.model import Cohort, Patient, Status, VitalStatus, Disease, Age
from gpsea.io import GpseaJSONDecoder
from gpsea.analysis import StatisticResult
from gpsea.analysis.predicate.genotype import (
    GenotypePolyPredicate,
    VariantPredicates,
    diagnosis_predicate,
    monoallelic_predicate,
)
from gpsea.analysis.temporal import SurvivalAnalysis, SurvivalAnalysisResult, Survival
from gpsea.analysis.temporal.endpoint import hpo_onset, death
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
        a_predicate=in_exon_3,
        b_predicate=~in_exon_3,
        a_label="Exon 3",
        b_label="Other exon",
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

    def test_raises_an_exception_for_an_invalid_dataset(
        self,
        survival_analysis: SurvivalAnalysis,
    ):
        d_one = Disease.from_raw_parts("OMIM:100000", name="One", is_observed=True)
        d_two = Disease.from_raw_parts("OMIM:200000", name="Two", is_observed=True)
        gt_predicate = diagnosis_predicate(
            diagnoses=(d.identifier for d in (d_one, d_two)),
            labels=(d.name for d in (d_one, d_two)),
        )
        endpoint = death()

        died_at_1y = VitalStatus(
            status=Status.DECEASED,
            age_of_death=Age.from_iso8601_period("P1Y"),
        )
        cohort = (
            Patient.from_raw_parts("A", vital_status=died_at_1y, diseases=(d_one,)),
            Patient.from_raw_parts("B", vital_status=died_at_1y, diseases=(d_one,)),
            Patient.from_raw_parts("C", vital_status=died_at_1y, diseases=(d_one,)),
            Patient.from_raw_parts("D", vital_status=died_at_1y, diseases=(d_two,)),
            Patient.from_raw_parts("E", vital_status=died_at_1y, diseases=(d_two,)),
        )

        with pytest.raises(AnalysisException) as e:
            survival_analysis.compare_genotype_vs_survival(
                cohort=cohort,
                gt_predicate=gt_predicate,
                endpoint=endpoint,
            )

        assert e.value.args == (
            "The survival values did not meet the expectation of the statistical test!",
        )


class TestSurvivalAnalysisResult:

    @pytest.fixture(scope="class")
    def result(
        self,
        umod_gt_predicate: GenotypePolyPredicate,
    ) -> SurvivalAnalysisResult:
        data = pd.DataFrame(
            data={
                "patient_id": ["A", "B", "C", "D", "E", "F", "G", "H"],
                "genotype": [0, 1, None, 0, 0, 1, 1, 1],
                "phenotype": [
                    Survival(value=12.3, is_censored=False),
                    None,
                    Survival(value=45.6, is_censored=True),
                    Survival(value=22.6, is_censored=True),
                    Survival(value=45.6, is_censored=False),
                    Survival(value=30.0, is_censored=False),
                    Survival(value=68.0, is_censored=False),
                    Survival(value=95.0, is_censored=True),
                ],
            }
        ).set_index("patient_id")
        return SurvivalAnalysisResult(
            gt_predicate=umod_gt_predicate,
            endpoint=death(),
            statistic=LogRankTest(),
            data=data,
            statistic_result=StatisticResult(statistic=1.0, pval=0.1234),
        )

    def test_properties(
        self,
        result: SurvivalAnalysisResult,
    ):
        assert tuple(result.gt_predicate.group_labels) == (
            "Exon 3",
            "Other exon",
        )
        assert result.data.index.to_list() == ["A", "B", "C", "D", "E", "F", "G", "H"]
        assert result.pval == pytest.approx(0.1234)

    def test_complete_records(
        self,
        result: SurvivalAnalysisResult,
    ):
        records = result.complete_records()

        assert records.shape == (6, 2)
        assert records.loc["A", "genotype"] == 0
        assert records.loc["A", "phenotype"] == Survival(value=12.3, is_censored=False)

    @pytest.mark.skip("Only for manual debugging")
    def test_plot_kaplan_meier(
        self,
        result: SurvivalAnalysisResult,
    ):
        fig, ax = plt.subplots(figsize=(6, 4), dpi=120)
        result.plot_kaplan_meier_curves(
            ax=ax,
        )
        fig.savefig("example-kaplan-meier.png")
