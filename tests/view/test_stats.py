import hpotk
import math
import pytest

import pandas as pd

from gpsea.analysis.pcats import HpoTermAnalysisResult
from gpsea.analysis.pcats.stats import FisherExactTest
from gpsea.analysis.predicate.genotype import GenotypePolyPredicate
from gpsea.analysis.predicate.phenotype import HpoPredicate
from gpsea.analysis.mtc_filter import PhenotypeMtcResult
from gpsea.view import MtcStatsViewer


class TestStatsViewable:

    @pytest.fixture(scope='class')
    def hpo_term_analysis_result(
        self,
        hpo: hpotk.MinimalOntology,
        suox_gt_predicate: GenotypePolyPredicate,
    ) -> HpoTermAnalysisResult:
        is_arachnodactyly = HpoPredicate(
            hpo=hpo,
            query=hpotk.TermId.from_curie('HP:0001166'),  # Arachnodactyly
        )
        is_seizure = HpoPredicate(
            hpo=hpo,
            query=hpotk.TermId.from_curie('HP:0001250'),  # Seizure
        )
        return HpoTermAnalysisResult(
            gt_predicate=suox_gt_predicate,
            statistic=FisherExactTest(),
            mtc_correction='fdr_bh',
            pheno_predicates=(
                is_arachnodactyly,
                is_seizure,
            ),
            n_usable=(40, 20),
            all_counts=(
                pd.DataFrame(
                    data=[[10, 5], [10, 15]],
                    index=pd.Index(is_arachnodactyly.get_categories()),
                    columns=pd.Index(suox_gt_predicate.get_categories())
                ),
                pd.DataFrame(
                    data=[[5, 0], [5, 10]],
                    index=pd.Index(is_seizure.get_categories()),
                    columns=pd.Index(suox_gt_predicate.get_categories())
                ),
            ),
            pvals=(math.nan, .005,),
            corrected_pvals=(math.nan, .01),
            mtc_filter_name='Random MTC filter',
            mtc_filter_results=(
                PhenotypeMtcResult.fail("RMF01", "Not too interesting"),
                PhenotypeMtcResult.ok(),
            ),
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
            report.write(fh)
