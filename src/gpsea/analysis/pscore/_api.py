import abc
import math
import typing

import numpy as np
import pandas as pd

from gpsea.model import Patient
from gpsea.config import PALETTE_DATA, PALETTE_SPECIAL
from ..clf import GenotypeClassifier
from .stats import PhenotypeScoreStatistic

from .._base import MonoPhenotypeAnalysisResult, Statistic, StatisticResult
from .._partition import ContinuousPartitioning


class PhenotypeScorer(ContinuousPartitioning, metaclass=abc.ABCMeta):
    """
    `PhenotypeScorer` assigns the patient with a phenotype score.

    The score can be `math.nan` if it is not possible to compute the score for a patient.

    The scorer can be created by wrapping a scoring function (see :func:`~PhenotypeScorer.wrap_scoring_function`).
    """

    @staticmethod
    def wrap_scoring_function(
        func: typing.Callable[[Patient], float],
        name: str = "Custom Scoring Function",
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
        return FunctionPhenotypeScorer(
            name=name,
            func=func,
        )

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

    @property
    def name(self) -> str:
        return self._name

    @property
    def description(self) -> str:
        return "A custom function to compute the phenotype score"

    @property
    def variable_name(self) -> str:
        return "Phenotype score"

    def __init__(
        self,
        name: str,
        func: typing.Callable[[Patient], float],
    ):
        assert isinstance(name, str)
        self._name = name
        self._func = func

    def score(self, patient: Patient) -> float:
        """
        Apply the function to compute the phenotype score.
        """
        return self._func(patient)


class PhenotypeScoreAnalysisResult(MonoPhenotypeAnalysisResult):
    """
    `PhenotypeScoreAnalysisResult` is a container for :class:`PhenotypeScoreAnalysis` results.

    The :attr:`~gpsea.analysis.MonoPhenotypeAnalysisResult.data`
    property provides a data frame with phenotype score for each tested individual:

    ==========  ========  =========
    patient_id  genotype  phenotype
    ==========  ========  =========
    patient_1   0         1
    patient_2   0         nan
    patient_3   None      2
    patient_4   1         2
    ...         ...       ...
    ==========  ========  =========

    The DataFrame index includes the identifiers of the tested individuals and the values are stored
    in `genotype` and `phenotype` columns.

    The `genotype` includes the genotype category ID (:attr:`~gpsea.analysis.clf.PatientCategory.cat_id`)
    or `None` if the patient cannot be assigned into any genotype category.

    The `phenotype` contains a `float` with the phenotype score. A `NaN` value is used
    if the phenotype score is impossible to compute.
    """

    def __init__(
        self,
        gt_clf: GenotypeClassifier,
        phenotype: PhenotypeScorer,
        statistic: Statistic,
        data: pd.DataFrame,
        statistic_result: StatisticResult,
    ):
        super().__init__(gt_clf, phenotype, statistic, data, statistic_result)
        assert isinstance(phenotype, PhenotypeScorer)

        # Check that the provided genotype predicate defines the same categories
        # as those found in `data.`
        actual = set(
            int(val)
            for val in data[MonoPhenotypeAnalysisResult.GT_COL].unique()
            if val is not None and not math.isnan(val)
        )
        expected = set(c.cat_id for c in self._gt_clf.get_categories())
        assert actual == expected, "Mismatch in the genotype classes"

    def phenotype_scorer(self) -> PhenotypeScorer:
        """
        Get the scorer that computed the phenotype score.
        """
        # We are sure that `self._phenotype` is a `PhenotypeScorer`
        # because of the instance check in `__init__` and `PhenotypeScorer`
        # being a subclass of `Partitioning`.
        return self._phenotype  # type: ignore

    def _make_data_df(
        self,
    ) -> pd.DataFrame:
        # skip the patients with unassigned genotype group
        not_na = self._data.notna()
        not_na_gts = not_na.all(axis="columns")
        return self._data.loc[not_na_gts]

    def _make_x_and_tick_labels(
        self,
        data: pd.DataFrame,
    ) -> typing.Tuple[
        typing.Sequence[typing.Sequence[float]],
        typing.Sequence[str],
    ]:
        x = [
            data.loc[
                data[MonoPhenotypeAnalysisResult.GT_COL] == c.category.cat_id,
                MonoPhenotypeAnalysisResult.PH_COL,
            ].to_list()
            for c in self._gt_clf.get_categorizations()
        ]

        gt_cat_names = [c.category.name for c in self._gt_clf.get_categorizations()]

        return x, gt_cat_names

    def plot_boxplots(
        self,
        ax,
        colors: typing.Sequence[str] = PALETTE_DATA,
        median_color: str = PALETTE_SPECIAL,
        **boxplot_kwargs,
    ):
        """
        Draw box plot with distributions of phenotype scores for the genotype groups.

        :param ax: the Matplotlib :class:`~matplotlib.axes.Axes` to draw on.
        :param colors: a sequence with color palette for the box plot patches.
        :param median_color: a `str` with the color for the boxplot median line.
        :param boxplot_kwargs: arguments to pass into :func:`matplotlib.axes.Axes.boxplot` function.
        """
        data = self._make_data_df()

        x, gt_cat_names = self._make_x_and_tick_labels(data)
        patch_artist = boxplot_kwargs.pop("patch_artist", True)
        tick_labels = boxplot_kwargs.pop("tick_labels", gt_cat_names)

        bplot = ax.boxplot(
            x=x,
            patch_artist=patch_artist,
            tick_labels=tick_labels,
            **boxplot_kwargs,
        )

        # Set face colors of the boxes
        col_idxs = self._choose_palette_idxs(
            n_categories=self._gt_clf.n_categorizations(), n_colors=len(colors)
        )
        for patch, col_idx in zip(bplot["boxes"], col_idxs):
            patch.set_facecolor(colors[col_idx])

        for median in bplot["medians"]:
            median.set_color(median_color)

    def plot_violins(
        self,
        ax,
        colors: typing.Sequence[str] = PALETTE_DATA,
        **violinplot_kwargs,
    ):
        """
        Draw a violin plot with distributions of phenotype scores for the genotype groups.

        :param ax: the Matplotlib :class:`~matplotlib.axes.Axes` to draw on.
        :param colors: a sequence with color palette for the violin patches.
        :param violinplot_kwargs: arguments to pass into :func:`matplotlib.axes.Axes.violinplot` function.
        """
        data = self._make_data_df()

        x, gt_cat_names = self._make_x_and_tick_labels(data)

        showmeans = violinplot_kwargs.pop("showmeans", False)
        showextrema = violinplot_kwargs.pop("showextrema", False)

        parts = ax.violinplot(
            dataset=x,
            showmeans=showmeans,
            showextrema=showextrema,
            **violinplot_kwargs,
        )

        # quartile1, medians, quartile3 = np.percentile(x, [25, 50, 75], axis=1)
        quartile1 = [np.percentile(v, 25) for v in x]
        medians = [np.median(v) for v in x]
        quartile3 = [np.percentile(v, 75) for v in x]
        x = [sorted(val) for val in x]
        whiskers = np.array(
            [
                PhenotypeScoreAnalysisResult._adjacent_values(sorted_array, q1, q3)
                for sorted_array, q1, q3 in zip(x, quartile1, quartile3)
            ]
        )
        whiskers_min, whiskers_max = whiskers[:, 0], whiskers[:, 1]

        inds = np.arange(1, len(medians) + 1)
        ax.scatter(inds, medians, marker="o", color="white", s=30, zorder=3)
        ax.vlines(inds, quartile1, quartile3, color="k", linestyle="-", lw=5)
        ax.vlines(inds, whiskers_min, whiskers_max, color="k", linestyle="-", lw=1)

        ax.xaxis.set(
            ticks=np.arange(1, len(gt_cat_names) + 1),
            ticklabels=gt_cat_names,
        )

        col_idxs = self._choose_palette_idxs(
            n_categories=self._gt_clf.n_categorizations(), n_colors=len(colors)
        )
        for pc, color_idx in zip(parts["bodies"], col_idxs):
            pc.set(
                facecolor=colors[color_idx],
                edgecolor=None,
                alpha=1,
            )

    @staticmethod
    def _adjacent_values(vals, q1, q3):
        upper_adjacent_value = q3 + (q3 - q1) * 1.5
        upper_adjacent_value = np.clip(upper_adjacent_value, q3, vals[-1])

        lower_adjacent_value = q1 - (q3 - q1) * 1.5
        lower_adjacent_value = np.clip(lower_adjacent_value, vals[0], q1)
        return lower_adjacent_value, upper_adjacent_value

    def __eq__(self, value: object) -> bool:
        return isinstance(value, PhenotypeScoreAnalysisResult) and super(
            MonoPhenotypeAnalysisResult, self
        ).__eq__(value)

    def __hash__(self) -> int:
        return super(MonoPhenotypeAnalysisResult, self).__hash__()

    def __str__(self) -> str:
        return (
            "PhenotypeScoreAnalysisResult("
            f"gt_clf={self._gt_clf}, "
            f"phenotype_scorer={self._phenotype}, "
            f"statistic={self._statistic}, "
            f"data={self._data}, "
            f"statistic_result={self._statistic_result})"
        )

    def __repr__(self) -> str:
        return str(self)


class PhenotypeScoreAnalysis:
    """
    `PhenotypeScoreAnalysis` tests the association between two or more genotype classes
    and a phenotype score.

    A genotype class is assigned by a :class:`~gpsea.analysis.clf.GenotypeClassifier`
    and the phenotype score is computed with a :class:`~gpsea.analysis.pscore.PhenotypeScorer`.

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
        gt_clf: GenotypeClassifier,
        pheno_scorer: PhenotypeScorer,
    ) -> PhenotypeScoreAnalysisResult:
        """
        Compute the association between genotype groups and phenotype score.

        :param cohort: the cohort to analyze.
        :param gt_clf: a classifier for assigning an individual into a genotype class.
        :param pheno_scorer: the scorer to compute phenotype score.
        """
        assert (
            gt_clf.n_categorizations() == 2
        ), "We only support 2 genotype categories at this point"
        assert isinstance(pheno_scorer, PhenotypeScorer)

        idx = pd.Index((patient.patient_id for patient in cohort), name="patient_id")
        data = pd.DataFrame(
            None,
            index=idx,
            columns=MonoPhenotypeAnalysisResult.DATA_COLUMNS,
        )

        # Apply the classifier and scorer on the individuals
        for individual in cohort:
            gt_cat = gt_clf.test(individual)
            if gt_cat is None:
                data.loc[
                    individual.patient_id,
                    MonoPhenotypeAnalysisResult.GT_COL
                ] = None
            else:
                data.loc[
                    individual.patient_id, 
                    MonoPhenotypeAnalysisResult.GT_COL
                ] = gt_cat.category.cat_id

            data.loc[
                individual.patient_id,
                MonoPhenotypeAnalysisResult.PH_COL
            ] = pheno_scorer.score(individual)

        # Sort by PatientCategory.cat_id and unpack.
        # For now, we only allow to have up to 2 groups.
        x_key, y_key = sorted(
            data[MonoPhenotypeAnalysisResult.GT_COL].dropna().unique()
        )
        x = data.loc[
            data[MonoPhenotypeAnalysisResult.GT_COL] == x_key,
            MonoPhenotypeAnalysisResult.PH_COL,
        ].to_numpy(dtype=float)  # type: ignore
        y = data.loc[
            data[MonoPhenotypeAnalysisResult.GT_COL] == y_key,
            MonoPhenotypeAnalysisResult.PH_COL,
        ].to_numpy(dtype=float)  # type: ignore
        result = self._statistic.compute_pval(scores=(x, y))

        return PhenotypeScoreAnalysisResult(
            gt_clf=gt_clf,
            phenotype=pheno_scorer,
            statistic=self._statistic,
            data=data,
            statistic_result=result,
        )
