import pytest

from gpsea.analysis import HpoMtcReport
from gpsea.view import MtcStatsViewer


class TestStatsViewable:

    @pytest.fixture
    def stats_viewer(self) -> MtcStatsViewer:
        return MtcStatsViewer()

    @pytest.mark.skip('Until we design a more reasonable test')
    def test_process(
        self,
        stats_viewer: MtcStatsViewer,
    ):
        mtc_report = HpoMtcReport(
            filter_name='identity filter',
            mtc_name='bonferroni',
            filter_results_map={
                # The reason for skipping a phenotype -> the number of phenotypes skipped for the reason
                'I slept bad tonight': 0,
                'The sun is too bright today': 5,
                'Life is a conspiracy': 80,
                'I need coffee': 7,
            },
            n_terms_before_filtering=100,  # The filtered out (80 + 7 + 5) + the unfiltered
        )

        report = stats_viewer.process(hpo_mtc_report=mtc_report)
        with open('mtc_stats.html', 'w') as fh:
            fh.write(report)
