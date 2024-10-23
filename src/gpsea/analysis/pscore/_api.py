import abc
import typing

import pandas as pd

from gpsea.model import Patient
from ..predicate.genotype import GenotypePolyPredicate
from .stats import PhenotypeScoreStatistic

from .._base import MonoPhenotypeAnalysisResult


class PhenotypeScorer(metaclass=abc.ABCMeta):
    """
    `PhenotypeScorer` assigns the patient with a phenotype score.

    The score can be `math.nan` if it is not possible to compute the score for a patient.

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

    @abc.abstractmethod
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


class PhenotypeScoreAnalysisResult(MonoPhenotypeAnalysisResult):
    """
    `PhenotypeScoreAnalysisResult` is a container for :class:`PhenotypeScoreAnalysis` results.

    The :attr:`PhenotypeScoreAnalysisResult.data` property provides a data frame
    with phenotype score for each tested individual:
    
    ==========  ========  =========
    patient_id  genotype  phenotype
    ==========  ========  =========
    patient_1   0         1
    patient_2   0         nan
    patient_3   None      2
    patient_4   1         2
    ...         ...       ...
    ==========  ========  ===============

    The DataFrame index includes the identifiers of the tested individuals and the values are stored
    in `genotype` and `phenotype` columns.

    The `genotype` includes the genotype category ID (:attr:`~gpsea.analysis.predicate.PatientCategory.cat_id`)
    or `None` if the patient cannot be assigned into any genotype category.
    
    The `phenotype` contains a `float` with the phenotype score. A `NaN` value is used
    if the phenotype score is impossible to compute.
    """

    def plot_boxplots(
        self,
        ax,
        colors=("darksalmon", "honeydew"),
    ):
        """
        Draw box plot with distributions of phenotype scores for the genotype groups.

        :param gt_predicate: the genotype predicate used to produce the genotype groups.
        :param ax: the Matplotlib :class:`~matplotlib.axes.Axes` to draw on.
        :param colors: a tuple with colors to use for coloring the box patches of the box plot.
        """
        # skip the patients with unassigned genotype group
        bla = self._data.notna()
        not_na_gts = bla.all(axis='columns')
        data = self._data.loc[not_na_gts]
        
        # Check that the provided genotype predicate defines the same categories
        # as those found in `data.`
        actual = set(data[MonoPhenotypeAnalysisResult.GT_COL].unique())
        expected = set(c.cat_id for c in self._gt_predicate.get_categories())
        assert actual == expected, 'Mismatch in the genotype categories'
        
        x = [
            data.loc[
                data[MonoPhenotypeAnalysisResult.GT_COL] == c.category.cat_id, MonoPhenotypeAnalysisResult.PH_COL
            ].to_list()
            for c in self._gt_predicate.get_categorizations()
        ]
        
        gt_cat_names = [
            c.category.name for c in self._gt_predicate.get_categorizations()
        ]
        bplot = ax.boxplot(
            x=x,
            patch_artist=True,
            tick_labels=gt_cat_names,
        )

        for patch, color in zip(bplot["boxes"], colors):
            patch.set_facecolor(color)

    def __eq__(self, value: object) -> bool:
        return isinstance(value, PhenotypeScoreAnalysisResult) \
            and super(MonoPhenotypeAnalysisResult, self).__eq__(value)
    
    def __hash__(self) -> int:
        return super(MonoPhenotypeAnalysisResult, self).__hash__()

    def __str__(self) -> str:
        return (
            "PhenotypeScoreAnalysisResult("
            "gt_predicate={self._gt_predicate}, "
            "statistic={self._statistic}, "
            "data={self._data}, "
            "pval={self._pval})"
        )

    def __repr__(self) -> str:
        return str(self)


class PhenotypeScoreAnalysis:
    """
    `PhenotypeScoreAnalysis` tests the association between two or more genotype groups
    and a phenotype score.

    The genotype groups are created by a :class:`~gpsea.analysis.predicate.genotype.GenotypePolyPredicate`
    and the phenotype score is computed with :class:`~gpsea.analysis.pscore.PhenotypeScorer`.

    The association is tested with a :class:`~gpsea.analysis.pscore.stats.PhenotypeScoreStatistic`
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
            columns=MonoPhenotypeAnalysisResult.DATA_COLUMNS,
        )

        # Apply the predicates on the patients
        for patient in cohort:
            gt_cat = gt_predicate.test(patient)
            if gt_cat is None:
                data.loc[patient.patient_id, MonoPhenotypeAnalysisResult.GT_COL] = None
            else:
                data.loc[patient.patient_id, MonoPhenotypeAnalysisResult.GT_COL] = gt_cat.category.cat_id
            
            data.loc[patient.patient_id, MonoPhenotypeAnalysisResult.PH_COL] = pheno_scorer.score(patient)

        # Sort by PatientCategory.cat_id and unpack.
        # For now, we only allow to have up to 2 groups.
        x_key, y_key = sorted(data[MonoPhenotypeAnalysisResult.GT_COL].dropna().unique())
        x = data.loc[
            data[MonoPhenotypeAnalysisResult.GT_COL] == x_key, MonoPhenotypeAnalysisResult.PH_COL
        ].to_numpy(dtype=float)  # type: ignore
        y = data.loc[
            data[MonoPhenotypeAnalysisResult.GT_COL] == y_key, MonoPhenotypeAnalysisResult.PH_COL
        ].to_numpy(dtype=float)  # type: ignore
        pval = self._statistic.compute_pval(scores=(x, y))

        return PhenotypeScoreAnalysisResult(
            gt_predicate=gt_predicate,
            statistic=self._statistic,
            data=data,
            pval=pval,
        )
