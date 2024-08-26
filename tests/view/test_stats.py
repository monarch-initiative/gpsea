import hpotk
import math
import pytest

import pandas as pd

from gpsea.analysis.multip import HpoTermAnalysisResult
from gpsea.analysis.predicate import PatientCategories
from gpsea.analysis.predicate.genotype import GenotypePolyPredicate
from gpsea.analysis.mtc_filter import PhenotypeMtcResult
from gpsea.view import MtcStatsViewer


class TestStatsViewable:

    @pytest.fixture(scope='class')
    def hpo_term_analysis_result(
        self,
        suox_gt_predicate: GenotypePolyPredicate,
    ) -> HpoTermAnalysisResult:
        return HpoTermAnalysisResult(
            phenotypes=(
                hpotk.TermId.from_curie('HP:0001166'),  # Arachnodactyly
                hpotk.TermId.from_curie('HP:0001250'),  # Seizure
            ),
            n_usable=(40, 20),
            all_counts=(
                pd.DataFrame(
                    data=[[10, 5], [10, 15]],
                    index=pd.Index((PatientCategories.YES, PatientCategories.NO,)),
                    columns=pd.Index(suox_gt_predicate.get_categories())
                ),
                pd.DataFrame(
                    data=[[5, 0], [5, 10]],
                    index=pd.Index((PatientCategories.YES, PatientCategories.NO,)),
                    columns=pd.Index(suox_gt_predicate.get_categories())
                ),
            ),
            pvals=(math.nan, .005,),
            corrected_pvals=(math.nan, .01),
            gt_predicate=suox_gt_predicate,
            mtc_filter_name='Random MTC filter',
            mtc_filter_results=(
                PhenotypeMtcResult.fail("Not too interesting"),
                PhenotypeMtcResult.ok(),
            ),
            mtc_name='fdr_bh',
        )

    @pytest.fixture
    def stats_viewer(self) -> MtcStatsViewer:
        return MtcStatsViewer()

    @pytest.mark.skip('Until we design a more reasonable test')
    def test_process(
        self,
        stats_viewer: MtcStatsViewer,
        hpo_term_analysis_result: HpoTermAnalysisResult
    ):
        report = stats_viewer.process(result=hpo_term_analysis_result)
        with open('mtc_stats.html', 'w') as fh:
            fh.write(report)
