import math

import hpotk
import pandas as pd
import pytest

from gpsea.analysis import MultiPhenotypeAnalysisResult, StatisticResult
from gpsea.analysis.predicate.genotype import GenotypePolyPredicate
from gpsea.analysis.predicate.phenotype import HpoPredicate
from gpsea.analysis.pcats.stats import FisherExactTest


@pytest.fixture(scope="class")
def multi_phenotype_analysis_result(
    hpo: hpotk.MinimalOntology,
    suox_gt_predicate: GenotypePolyPredicate,
) -> MultiPhenotypeAnalysisResult:
    is_arachnodactyly = HpoPredicate(
        hpo=hpo,
        query=hpotk.TermId.from_curie("HP:0001166"),  # Arachnodactyly
    )
    is_seizure = HpoPredicate(
        hpo=hpo,
        query=hpotk.TermId.from_curie("HP:0001250"),  # Seizure
    )
    is_polydactyly = HpoPredicate(
        hpo=hpo,
        query=hpotk.TermId.from_curie("HP:0100259"),  # Postaxial polydactyly
    )
    is_clinodactyly = HpoPredicate(
        hpo=hpo,
        query=hpotk.TermId.from_curie("HP:0030084"),  # Clinodactyly
    )
    return MultiPhenotypeAnalysisResult(
        gt_predicate=suox_gt_predicate,
        statistic=FisherExactTest(),
        mtc_correction="fdr_bh",
        pheno_predicates=(
            is_arachnodactyly,
            is_seizure,
            is_polydactyly,
            is_clinodactyly,
        ),
        n_usable=(40, 20, 115, 10),
        all_counts=(
            pd.DataFrame(
                data=[[10, 5], [10, 15]],
                index=pd.Index(is_arachnodactyly.get_categories()),
                columns=pd.Index(suox_gt_predicate.get_categories()),
            ),
            pd.DataFrame(
                data=[[5, 0], [5, 10]],
                index=pd.Index(is_seizure.get_categories()),
                columns=pd.Index(suox_gt_predicate.get_categories()),
            ),
            pd.DataFrame(
                data=[[50, 0], [5, 60]],
                index=pd.Index(is_polydactyly.get_categories()),
                columns=pd.Index(suox_gt_predicate.get_categories()),
            ),
            pd.DataFrame(
                data=[[0, 0], [10, 0]],
                index=pd.Index(is_clinodactyly.get_categories()),
                columns=pd.Index(suox_gt_predicate.get_categories()),
            ),
        ),
        statistic_results=(
            None,
            StatisticResult(statistic=1., pval=0.005),
            StatisticResult(statistic=10.0, pval=0.0005),
            StatisticResult(statistic=0.1, pval=0.05),
        ),
        corrected_pvals=(math.nan, 0.01, 0.001, 0.5),
    )


class TestMultiPhenotypeAnalysisResult:

    def test_significant_phenotype_indices(
        self,
        multi_phenotype_analysis_result: MultiPhenotypeAnalysisResult,
    ):
        indices = multi_phenotype_analysis_result.significant_phenotype_indices(
            alpha=0.05, pval_kind="corrected"
        )
        assert indices == (2, 1)
