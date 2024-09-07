import hpotk
import pytest

from gpsea.analysis.pcats import HpoTermAnalysisResult
from gpsea.view import HpoTermAnalysisResultFormatter


class TestHpoTermAnalysisResultFormatter:

    @pytest.fixture
    def formatter(
        self,
        hpo: hpotk.MinimalOntology,
    ) -> HpoTermAnalysisResultFormatter:
        return HpoTermAnalysisResultFormatter(
            hpo=hpo,
        )

    @pytest.mark.skip("Only for manual control")
    def test_summarize(
        self,
        formatter: HpoTermAnalysisResultFormatter,
        hpo_result: HpoTermAnalysisResult,
    ):
        df = formatter.make_summary_dataframe(result=hpo_result)

        print(df)
