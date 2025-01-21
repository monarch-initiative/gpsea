import hpotk
import pytest

from gpsea.analysis.pcats import HpoTermAnalysisResult
from gpsea.view import summarize_hpo_analysis


@pytest.mark.skip("Just for manual testing and debugging")
def test_summarize_hpo_analysis(
    hpo: hpotk.MinimalOntology,
    hpo_term_analysis_result: HpoTermAnalysisResult,
):
    df = summarize_hpo_analysis(
        hpo=hpo,
        result=hpo_term_analysis_result,
    )

    print(df)
