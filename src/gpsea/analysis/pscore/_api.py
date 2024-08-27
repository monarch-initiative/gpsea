import abc
import typing

import pandas as pd

from gpsea.model import Patient
from ..predicate.genotype import GenotypePolyPredicate
from .stats import PhenotypeScoreStatistic


class PhenotypeScorer(metaclass=abc.ABCMeta):
    """
    `PhenotypeScorer` assigns the patient with a phenotype score.

    The score can be :attr:`math.nan` if it is not possible to compute the score for a patient.

    The scorer can be created by wrapping a scoring function (see :func:`~PhenotypeScorer.wrap_scoring_function`).
    """

    @staticmethod
    def wrap_scoring_function(
        func: typing.Callable[[Patient], float],
    ) -> "PhenotypeScorer":
        """
        Create a `PhenotypeScorer` by wrap the provided scoring function `func`.

        The function must take exactly one argument of type :class:`~gpsea.model.Patient`
        and return a `float` with the corresponding phenotype score.

        Example
        ^^^^^^^

        >>> from gpsea.analysis.pscore import PhenotypeScorer
        >>> def f(p): 123.4
        >>> phenotype_scorer = PhenotypeScorer.wrap_scoring_function(f)

        `phenotype_scorer` will assign all patients a score of `123.4`.

        :param func: the scoring function.
        """
        return FunctionPhenotypeScorer(func=func)

    def score(self, patient: Patient) -> float:
        """
        Compute the score for the `patient`.
        """
        pass


class FunctionPhenotypeScorer(PhenotypeScorer):
    """
    `FunctionPhenotypeScorer` computes the phenotype score using the provided function/closure.
    """
    # NOT PART OF THE PUBLIC API

    def __init__(
        self,
        func: typing.Callable[[Patient], float],
    ):
        self._func = func

    def score(self, patient: Patient) -> float:
        """
        Apply the function to compute the phenotype score.
        """
        return self._func(patient)


class PhenotypeScoreAnalysisResult:
    """
    `PhenotypeScoreAnalysisResult` is a container for :class:`PhenotypeScoreAnalysis` results.
    """

    def __init__(
        self,
        genotype_phenotype_scores: pd.DataFrame,
        pval: float,
    ):
        self._genotype_phenotype_scores = genotype_phenotype_scores
        self._pval = float(pval)

    @property
    def genotype_phenotype_scores(self) -> pd.DataFrame:
        """
        Get the DataFrame with the genotype group and the phenotype score for each patient.

        The DataFrame has the following structure:

        ==========  ========  =========
        patient_id  genotype  phenotype
        ==========  ========  =========
        patient_1   0         1
        patient_2   0         3
        patient_3   None      2
        patient_4   1         2
        ...         ...       ...
        ==========  ========  =========

        The DataFrame index includes the patient IDs, and then there are 2 columns
        with the `genotype` group id (:attr:`~gpsea.analysis.predicate.PatientCategory.cat_id`)
        and the `phenotype` score. A `genotype` value may be missing if the patient
        cannot be assigned into any genotype category.
        """
        return self._genotype_phenotype_scores

    @property
    def pval(self) -> float:
        """
        Get the p value of the test.
        """
        return self._pval

    def plot_boxplots(
        self,
        gt_predicate: GenotypePolyPredicate,
        ax,
        colors=["darksalmon", "honeydew"],
    ):
        """
        Draw box plots with distributions of phenotype scores for genotype groups
        """
        # skip the patients with unassigned genotype group
        not_na_gts = self._genotype_phenotype_scores["genotype"].notna()
        data = self._genotype_phenotype_scores.loc[not_na_gts]
        
        # Check that the provided genotype predicate defines the same categories
        # as those found in `data.`
        actual = set(data["genotype"].unique())
        expected = set(c.cat_id for c in gt_predicate.get_categories())
        assert actual == expected, 'Mismatch in the genotype categories'
        
        x = [
            data.loc[data["genotype"] == c.category.cat_id, "phenotype"].to_list()
            for c in gt_predicate.get_categorizations()
        ]
        
        gt_cat_names = [
            c.category.name for c in gt_predicate.get_categorizations()
        ]
        bplot = ax.boxplot(
            x=x,
            patch_artist=True,
            tick_labels=gt_cat_names,
        )

        for patch, color in zip(bplot["boxes"], colors):
            patch.set_facecolor(color)


class PhenotypeScoreAnalysis:
    """
    `PhenotypeScoreAnalysis` tests the association between two or more genotype groups
    and a phenotype score.

    The genotype groups are created by a :class:`~gpsea.analysis.predicate.genotype.GenotypePolyPredicate`
    and the phenotype score is computed with :class:`~gpsea.analysis.pscore.PhenotypeScorer`.

    The association is tested with a :class:`~gpsea.analysis.pscore.PhenotypeScoreStatistic`
    and the results are reported as a :class:`PhenotypeScoreAnalysisResult`.
    """

    def __init__(
        self,
        score_statistic: PhenotypeScoreStatistic,
    ):
        assert isinstance(score_statistic, PhenotypeScoreStatistic)
        self._statistic = score_statistic

    def compare_genotype_vs_phenotype_score(
        self,
        cohort: typing.Iterable[Patient],
        gt_predicate: GenotypePolyPredicate,
        pheno_scorer: PhenotypeScorer,
    ) -> PhenotypeScoreAnalysisResult:
        """
        Compute the association between genotype groups and phenotype score.

        :param cohort: the cohort to analyze.
        :param gt_predicate: a predicate for assigning an individual into a genotype group.
        :param pheno_scorer: the scorer to compute phenotype score.
        """
        assert (
            gt_predicate.n_categorizations() == 2
        ), "We only support 2 genotype categories at this point"

        idx = pd.Index((patient.patient_id for patient in cohort), name="patient_id")
        data = pd.DataFrame(
            None,
            index=idx,
            columns=["genotype", "phenotype"],
        )

        # Apply the predicates on the patients
        for patient in cohort:
            gt_cat = gt_predicate.test(patient)
            if gt_cat is None:
                data.loc[patient.patient_id, "genotype"] = None
            else:
                data.loc[patient.patient_id, "genotype"] = gt_cat.category.cat_id
            
            data.loc[patient.patient_id, "phenotype"] = pheno_scorer.score(patient)

        # Sort by PatientCategory.cat_id and unpack.
        # For now, we only allow to have up to 2 groups.
        x_key, y_key = sorted(data["genotype"].dropna().unique())
        x = data.loc[data["genotype"] == x_key, "phenotype"].to_numpy(dtype=float)
        y = data.loc[data["genotype"] == y_key, "phenotype"].to_numpy(dtype=float)
        pval = self._statistic.compute_pval(scores=(x, y))

        return PhenotypeScoreAnalysisResult(
            genotype_phenotype_scores=data,
            pval=pval,
        )
