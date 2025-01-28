import math

import hpotk
import pytest
import pandas as pd

from gpsea.analysis import StatisticResult
from gpsea.analysis.clf import HpoClassifier, GenotypeClassifier
from gpsea.analysis.mtc_filter import PhenotypeMtcResult
from gpsea.analysis.pcats import HpoTermAnalysisResult
from gpsea.analysis.pcats.stats import FisherExactTest


@pytest.fixture(scope="package")
def hpo_term_analysis_result(
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
