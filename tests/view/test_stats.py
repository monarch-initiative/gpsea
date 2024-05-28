import pytest

from genophenocorr.analysis import HpoMtcReport
from genophenocorr.view import StatsViewable


class TestStatsViewable:

    @pytest.fixture
    def stats_viewable(self) -> StatsViewable:
        return StatsViewable()

    @pytest.mark.skip('Until we design a more reasonable test')
    def test_process(
            self,
            stats_viewable: StatsViewable,
    ):
        mtc_report = HpoMtcReport(
            filter_name='identity filter',
            mtc_name='bonferroni',
            filter_results_map={
                # The reason for skipping a phenotype -> the number of phenotypes skipped for the reason
                'Skipping test because I slept bad tonight': 0,
                'The sun is too bright today': 5,
                'The kids are too noisy': 80,
                'I need coffee': 7,
            },
            term_count=100,  # The filtered out (80 + 7 + 5) + the unfiltered
        )

        report = stats_viewable.process(hpo_mtc_report=mtc_report)
        with open('stats_viewable.html', 'w') as fh:
            fh.write(report)
