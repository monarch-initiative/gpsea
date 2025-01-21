import math
import os

import hpotk
import pandas as pd
import pytest

from gpsea.analysis import StatisticResult
from gpsea.analysis.pcats import HpoTermAnalysisResult
from gpsea.analysis.pcats.stats import FisherExactTest
from gpsea.analysis.clf import GenotypeClassifier, HpoClassifier
from gpsea.analysis.mtc_filter import PhenotypeMtcResult
from gpsea.model import Cohort
from gpsea.view import (
    CohortViewer,
    CohortVariantViewer,
    MtcStatsViewer,
    summarize_hpo_analysis,
)


@pytest.mark.skip("Just for manual testing and debugging")
class TestCohortViewer:
    @pytest.fixture
    def cohort_viewer(
        self,
        hpo: hpotk.MinimalOntology,
    ) -> CohortViewer:
        return CohortViewer(
            hpo=hpo,
        )

    def test_process_suox_cohort(
        self,
        cohort_viewer: CohortViewer,
        suox_cohort: Cohort,
        suox_mane_tx_id: str,
    ):
        report = cohort_viewer.process(
            cohort=suox_cohort,
            transcript_id=suox_mane_tx_id,
        )

        with open(os.path.join("dev", "SUOX.cohort.html"), "w") as fh:
            report.write(fh)

    def test_process_cyp21a2_cohort(
        self,
        cohort_viewer: CohortViewer,
        cyp21a2_cohort: Cohort,
        cyp21a2_mane_tx_id: str,
    ):
        report = cohort_viewer.process(
            cohort=cyp21a2_cohort,
            transcript_id=cyp21a2_mane_tx_id,
        )

        with open(os.path.join("dev", "CYP21A2.cohort.html"), "w") as fh:
            report.write(fh)


@pytest.mark.skip("Just for manual testing and debugging")
def test_viewer(
    suox_mane_tx_id: str,
    suox_cohort: Cohort,
):
    viewer = CohortVariantViewer(tx_id=suox_mane_tx_id)
    html = viewer.process(suox_cohort)

    with open("all_variants.html", "w") as fh:
        html.write(fh)


class TestMtcStatsViewer:
    @pytest.fixture(scope="class")
    def hpo_term_analysis_result(
        self,
        hpo: hpotk.MinimalOntology,
        suox_gt_clf: GenotypeClassifier,
    ) -> HpoTermAnalysisResult:
        is_arachnodactyly = HpoClassifier(
            hpo=hpo,
            query=hpotk.TermId.from_curie("HP:0001166"),  # Arachnodactyly
        )
        is_seizure = HpoClassifier(
            hpo=hpo,
            query=hpotk.TermId.from_curie("HP:0001250"),  # Seizure
        )
        return HpoTermAnalysisResult(
            gt_clf=suox_gt_clf,
            statistic=FisherExactTest(),
            mtc_correction="fdr_bh",
            pheno_clfs=(
                is_arachnodactyly,
                is_seizure,
            ),
            n_usable=(40, 20),
            all_counts=(
                pd.DataFrame(
                    data=[[10, 5], [10, 15]],
                    index=pd.Index(is_arachnodactyly.get_categories()),
                    columns=pd.Index(suox_gt_clf.get_categories()),
                ),
                pd.DataFrame(
                    data=[[5, 0], [5, 10]],
                    index=pd.Index(is_seizure.get_categories()),
                    columns=pd.Index(suox_gt_clf.get_categories()),
                ),
            ),
            statistic_results=(
                StatisticResult(statistic=None, pval=math.nan),
                StatisticResult(statistic=1.23, pval=0.01),
            ),
            corrected_pvals=(math.nan, 0.01),
            mtc_filter_name="Random MTC filter",
            mtc_filter_results=(
                PhenotypeMtcResult.fail("RMF01", "Not too interesting"),
                PhenotypeMtcResult.ok(),
            ),
        )

    @pytest.fixture
    def stats_viewer(self) -> MtcStatsViewer:
        return MtcStatsViewer()

    @pytest.mark.skip("Just for manual testing and debugging")
    def test_process(
        self,
        stats_viewer: MtcStatsViewer,
        hpo_term_analysis_result: HpoTermAnalysisResult,
    ):
        report = stats_viewer.process(result=hpo_term_analysis_result)
        with open("mtc_stats.html", "w") as fh:
            report.write(fh)


@pytest.mark.skip("Just for manual testing and debugging")
def test_summarize(
    hpo: hpotk.MinimalOntology,
    hpo_result: HpoTermAnalysisResult,
):
    df = summarize_hpo_analysis(hpo=hpo, result=hpo_result)

    print(df)
