import abc
import typing

from collections import defaultdict

import pandas as pd
import scipy.stats

from gpsea.model import Patient
from ..predicate.genotype import GenotypePolyPredicate

from ._base import Survival
from ._util import prepare_censored_data
from .stats import SurvivalStatistic

from .._base import MonoPhenotypeAnalysisResult
from .._partition import ContinuousPartitioning


class Endpoint(ContinuousPartitioning, metaclass=abc.ABCMeta):
    """
    `Endpoint` computes survival for the analyzed individual.

    An example endpoint includes :func:`~gpsea.analysis.survival.endpoint.death`,
    :func:`~gpsea.analysis.survival.endpoint.disease_onset`,
    or onset of a phenotypic feature (:func:`~gpsea.analysis.survival.endpoint.hpo_onset`).
    """

    @abc.abstractmethod
    def compute_survival(
        self,
        patient: Patient,
    ) -> typing.Optional[Survival]:
        """
        Compute a survival for a given `patient` or `None` if the `patient` lacks the required
        data (e.g. age of death or age at last investigation).
        """
        pass


class SurvivalAnalysisResult(MonoPhenotypeAnalysisResult):
    """
    `SurvivalAnalysisResult` includes the results of a :class:`~gpsea.analysis.temporal.SurvivalAnalysis`.

    The genotype categories and survival are reported in the `data` data frame with the following structure:

    ============  ===========  ============================================
     patient_id    genotype     phenotype
    ============  ===========  ============================================
     patient_1     0            `Survival(value=123.4, is_censored=False)`
     patient_2     0            `None`
     patient_3     `None`       `Survival(value=456.7, is_censored=True)`
     patient_4     1            `None`
     ...           ...          ...
    ============  ===========  ============================================

    The index includes the individual IDs (`patient_id`), and then there are 2 columns
    with the `genotype` group id (:attr:`~gpsea.analysis.predicate.PatientCategory.cat_id`)
    and the `phenotype` with the survival represented as :class:`~gpsea.analysis.tempo.Survival` object.

    A `genotype` value may be missing (`None`) if the individual cannot be assigned
    into a genotype category.
    Similarly, a `survival` is `None` if computing the survival for an individual is impossible.
    """

    def __init__(
        self,
        gt_predicate: GenotypePolyPredicate,
        endpoint: Endpoint,
        statistic: SurvivalStatistic,
        data: pd.DataFrame,
        pval: float,
    ):
        super().__init__(
            gt_predicate=gt_predicate,
            phenotype=endpoint,
            statistic=statistic,
            data=data,
            pval=pval,
        )
        assert isinstance(endpoint, Endpoint)

    @property
    def endpoint(self) -> Endpoint:
        """
        Get the endpoint used to compute the survival of the individuals.
        """
        # We are sure that `self._phenotype` is assignable to `Endpoint`
        # because of the instance check in `__init__` and `Endpoint`
        # being a subclass of `Partitioning`.
        return self._phenotype  # type: ignore

    def plot_kaplan_meier_curves(
        self,
        ax,
    ):
        """
        Plot genotype group survivals on the provided axes.

        The axes includes legend. However, if no survival is available
        for a genotype group, the group name will be missing from the legend.

        :param ax: a Matplotlib `Axes` to draw on.
        """
        for pat_cat in self._gt_predicate.get_categories():
            survivals = self._data.loc[
                self._data[MonoPhenotypeAnalysisResult.GT_COL] == pat_cat.cat_id,
                MonoPhenotypeAnalysisResult.PH_COL,
            ]
            non_na = survivals[survivals.notna()]
            if len(non_na) > 0:
                censored_data = prepare_censored_data(survivals=non_na)
                data = scipy.stats.ecdf(censored_data)
                data.sf.plot(ax, label=pat_cat.name)
        
        ax.legend()

    def __eq__(self, value: object) -> bool:
        return (
            isinstance(value, SurvivalAnalysisResult)
            and super(MonoPhenotypeAnalysisResult, self).__eq__(value)
        )

    def __hash__(self) -> int:
        return super(MonoPhenotypeAnalysisResult, self).__hash__()

    def __str__(self) -> str:
        return (
            "SurvivalAnalysisResult("
            f"gt_predicate={self._gt_predicate}, "
            f"endpoint={self._phenotype}, "
            f"statistic={self._statistic}, "
            f"data={self._data}, "
            f"pval={self._pval})"
        )

    def __repr__(self) -> str:
        return str(self)


class SurvivalAnalysis:
    """
    `SurvivalAnalysis` compares the survivals of genotype groups with respect
    to an :class:`~gpsea.analysis.temporal.Endpoint`.
    
    The cohort is partitioned into groups using a genotype predicate
    and survival is computed for each cohort member. The difference between
    survivals is tested with selected :class:`~gpsea.analysis.temporal.stats.SurvivalStatistic`.
    """

    def __init__(
        self,
        statistic: SurvivalStatistic,
    ):
        assert isinstance(statistic, SurvivalStatistic)
        self._statistic = statistic

    def compare_genotype_vs_survival(
        self,
        cohort: typing.Iterable[Patient],
        gt_predicate: GenotypePolyPredicate,
        endpoint: Endpoint,
    ) -> SurvivalAnalysisResult:
        """
        Execute the survival analysis on a given `cohort`.
        """
        
        idx = pd.Index((patient.patient_id for patient in cohort), name="patient_id")
        data = pd.DataFrame(
            None,
            index=idx,
            columns=MonoPhenotypeAnalysisResult.DATA_COLUMNS,
        )
        survivals = defaultdict(list)
        # Apply the predicate and the survival metric on the cohort
        for patient in cohort:
            gt_cat = gt_predicate.test(patient)
            if gt_cat is None:
                data.loc[patient.patient_id, MonoPhenotypeAnalysisResult.GT_COL] = None
            else:
                data.loc[patient.patient_id, MonoPhenotypeAnalysisResult.GT_COL] = gt_cat.category.cat_id

            survival = endpoint.compute_survival(patient)
            data.loc[patient.patient_id, MonoPhenotypeAnalysisResult.PH_COL] = survival  # type: ignore

            if gt_cat is not None and survival is not None:
                survivals[gt_cat].append(survival)

        vals = tuple(survivals[gt_cat] for gt_cat in gt_predicate.get_categorizations())
        pval = self._statistic.compute_pval(vals)

        return SurvivalAnalysisResult(
            gt_predicate=gt_predicate,
            endpoint=endpoint,
            statistic=self._statistic,
            data=data,
            pval=pval,
        )