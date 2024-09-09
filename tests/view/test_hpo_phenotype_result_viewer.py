import hpotk
import pytest

from gpsea.analysis.pcats import HpoTermAnalysisResult
from gpsea.view import HpoTermAnalysisResultViewer


class TestHpoTermAnalysisResultFormatter:

    @pytest.fixture
    def formatter(
        self,
        hpo: hpotk.MinimalOntology,
    ) -> HpoTermAnalysisResultViewer:
        return HpoTermAnalysisResultViewer(
            hpo=hpo,
        )

    @pytest.mark.skip("Only for manual control")
    def test_summarize(
        self,
        formatter: HpoTermAnalysisResultViewer,
        hpo_result: HpoTermAnalysisResult,
    ):
        df = formatter.make_summary_dataframe(result=hpo_result)

        print(df)
