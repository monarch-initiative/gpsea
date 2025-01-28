import math
import typing

import hpotk
import pandas as pd
import pytest

from gpsea.analysis import MultiPhenotypeAnalysisResult, StatisticResult
from gpsea.analysis.clf import GenotypeClassifier, HpoClassifier
from gpsea.analysis.pcats.stats import FisherExactTest


@pytest.fixture(scope="class")
def multi_phenotype_analysis_result(
    hpo: hpotk.MinimalOntology,
    suox_gt_clf: GenotypeClassifier,
) -> MultiPhenotypeAnalysisResult:
    is_arachnodactyly = HpoClassifier(
        hpo=hpo,
        query=hpotk.TermId.from_curie("HP:0001166"),  # Arachnodactyly
    )
    is_seizure = HpoClassifier(
        hpo=hpo,
        query=hpotk.TermId.from_curie("HP:0001250"),  # Seizure
    )
    is_polydactyly = HpoClassifier(
        hpo=hpo,
        query=hpotk.TermId.from_curie("HP:0100259"),  # Postaxial polydactyly
    )
    is_clinodactyly = HpoClassifier(
        hpo=hpo,
        query=hpotk.TermId.from_curie("HP:0030084"),  # Clinodactyly
    )
    return MultiPhenotypeAnalysisResult(
        gt_clf=suox_gt_clf,
        statistic=FisherExactTest(),
        mtc_correction="fdr_bh",
        pheno_clfs=(
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
                columns=pd.Index(suox_gt_clf.get_categories()),
            ),
            pd.DataFrame(
                data=[[5, 0], [5, 10]],
                index=pd.Index(is_seizure.get_categories()),
                columns=pd.Index(suox_gt_clf.get_categories()),
            ),
            pd.DataFrame(
                data=[[50, 0], [5, 60]],
                index=pd.Index(is_polydactyly.get_categories()),
                columns=pd.Index(suox_gt_clf.get_categories()),
            ),
            pd.DataFrame(
                data=[[0, 0], [10, 0]],
                index=pd.Index(is_clinodactyly.get_categories()),
                columns=pd.Index(suox_gt_clf.get_categories()),
            ),
        ),
        statistic_results=(
            None,
            StatisticResult(statistic=1.0, pval=0.005),
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

    @pytest.mark.parametrize(
        "n_cats, n_cols, expected",
        [
            [2, 2, (0, 1)],
            [2, 7, (0, 6)],
            [2, 10, (0, 9)],
            [2, 5, (0, 4)],
            [3, 5, (0, 2, 4)],
            [4, 5, (0, 1, 2, 4)],
            [5, 5, (0, 1, 2, 3, 4)],
        ],
    )
    def test__choose_palette_idxs(
        self,
        multi_phenotype_analysis_result: MultiPhenotypeAnalysisResult,
        n_cats: int,
        n_cols: int,
        expected: typing.Sequence[int],
    ):
        actual = multi_phenotype_analysis_result._choose_palette_idxs(
            n_categories=n_cats,
            n_colors=n_cols,
        )

        assert actual == expected

    def test__choose_palette_idxs__too_few_colors(
        self,
        multi_phenotype_analysis_result: MultiPhenotypeAnalysisResult,
    ):
        with pytest.raises(ValueError) as e:
            multi_phenotype_analysis_result._choose_palette_idxs(
                n_categories=2,
                n_colors=1,
            )
        assert e.value.args == ("Expected a palette with at least 2 colors but got 1",)

    def test__choose_palette_idxs__too_few_colors_for_cats(
        self,
        multi_phenotype_analysis_result: MultiPhenotypeAnalysisResult,
    ):
        with pytest.raises(ValueError) as e:
            multi_phenotype_analysis_result._choose_palette_idxs(
                n_categories=3,
                n_colors=2,
            )
        assert e.value.args == (
            "The predicate produces 3 categories but the palette includes only 2 colors!",
        )
