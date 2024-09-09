import hpotk
import pytest

from gpsea.analysis.pcats import HpoTermAnalysisResult
from gpsea.view import summarize_hpo_analysis


@pytest.mark.skip("Only for manual control")
def test_summarize(
    self,
    hpo: hpotk.MinimalOntology,
    hpo_result: HpoTermAnalysisResult,
):
    df = summarize_hpo_analysis(hpo=hpo, result=hpo_result)

    print(df)
