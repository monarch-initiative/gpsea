import pytest

from genophenocorr.analysis import HpoMtcReport
from genophenocorr.view import StatsViewer


class TestStatsViewable:

    @pytest.fixture
    def stats_viewer(self) -> StatsViewer:
        return StatsViewer()

    @pytest.mark.skip('Until we design a more reasonable test')
    def test_process(
            self,
            stats_viewer: StatsViewer,
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
            term_count=100,  # The filtered out (80 + 7 + 5) + the unfiltered
        )

        report = stats_viewer.process(hpo_mtc_report=mtc_report)
        with open('mtc_stats.html', 'w') as fh:
            fh.write(report)
