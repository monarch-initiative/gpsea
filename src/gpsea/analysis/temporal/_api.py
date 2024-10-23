import abc
import math
import typing

from collections import defaultdict

import pandas as pd

from gpsea.model import Patient
from ..predicate.genotype import GenotypePolyPredicate

from ._base import Survival
from .stats import SurvivalStatistic

from .._base import AnalysisResult


class Endpoint(metaclass=abc.ABCMeta):
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

    @abc.abstractmethod
    def display_question(self) -> str:
        pass


class SurvivalAnalysisResult(AnalysisResult):
    """
    `SurvivalAnalysisResult` includes the results of a :class:`~gpsea.analysis.temporal.SurvivalAnalysis`.
    """

    def __init__(
        self,
        gt_predicate: GenotypePolyPredicate,
        data: pd.DataFrame,
        pval: float,
    ):
        super().__init__(
            gt_predicate=gt_predicate,
        )

        assert isinstance(data, pd.DataFrame) and all(
            col in data for col in ("genotype", "survival")
        )
        self._data = data

        if isinstance(pval, float) and math.isfinite(pval) and 0.0 <= pval <= 1.0:
            self._pval = float(pval)
        else:
            raise ValueError(
                f"`p_val` must be a finite float in range [0, 1] but it was {pval}"
            )

    @property
    def data(self) -> pd.DataFrame:
        """
        Get the data frame with the genotype group
        and the corresponding :class:`~gpsea.analysis.tempo.Survival`.

        The DataFrame has the following structure:

        ==========  ==========  ============================================
        patient_id   genotype    survival
        ==========  ==========  ============================================
        patient_1   0            `Survival(value=123.4, is_censored=False)`
        patient_2   0            `None`
        patient_3   `None`       `Survival(value=456.7, is_censored=True)`
        patient_4   1            `None`
        ...         ...          ...
        ==========  ==========  ============================================

        The index includes the individual IDs (`patient_id`), and then there are 2 columns
        with the `genotype` group id (:attr:`~gpsea.analysis.predicate.PatientCategory.cat_id`)
        and the `survival` encoded as :class:`~gpsea.analysis.tempo.Survival` object.

        A `genotype` value may be missing (`None`) if the individual cannot be assigned
        into a genotype category.
        Similarly, a `survival` may be `None` if computing the survival is impossible for
        the individual in question.
        """
        return self._data

    def complete_records(self) -> pd.DataFrame:
        """
        Get the :meth:`~gpsea.analysis.temporal.SurvivalAnalysisResult.data` rows
        where both `genotype` and `survival` columns are available (i.e. not `None`).
        """
        return self._data.loc[
            self._data["genotype"].notna() & self._data["survival"].notna()
        ]

    @property
    def pval(self) -> float:
        """
        Get the p value of the test.
        """
        return self._pval

    def __eq__(self, value: object) -> bool:
        return (
            isinstance(value, SurvivalAnalysisResult)
            and super(AnalysisResult, self).__eq__(value)
            and self._data.equals(value._data)
            and self._pval == value._pval
        )

    def __hash__(self) -> int:
        return hash((
            super(AnalysisResult, self).__hash__(),
            self._data,
            self._pval,
        ))

    def __str__(self) -> str:
        return (
            "SurvivalAnalysisResult("
            "gt_predicate={self._gt_predicate}, "
            "data={self._data}, "
            "pval={self._pval})"
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
            columns=["genotype", "survival"],
        )
        survivals = defaultdict(list)
        # Apply the predicate and the survival metric on the cohort
        for patient in cohort:
            gt_cat = gt_predicate.test(patient)
            if gt_cat is None:
                data.loc[patient.patient_id, "genotype"] = None
            else:
                data.loc[patient.patient_id, "genotype"] = gt_cat.category.cat_id

            survival = endpoint.compute_survival(patient)
            data.loc[patient.patient_id, "survival"] = survival  # type: ignore

            if gt_cat is not None and survival is not None:
                survivals[gt_cat].append(survival)

        vals = tuple(survivals[gt_cat] for gt_cat in gt_predicate.get_categorizations())
        pval = self._statistic.compute_pval(vals)

        return SurvivalAnalysisResult(
            gt_predicate=gt_predicate,
            data=data,
            pval=pval,
        )
