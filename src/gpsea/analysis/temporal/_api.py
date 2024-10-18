import abc
import math
import typing

from collections import defaultdict

import pandas as pd

from gpsea.model import Patient
from ..predicate.genotype import GenotypePolyPredicate

from ._base import Survival
from .stats import SurvivalStatistic


class Endpoint(metaclass=abc.ABCMeta):

    @abc.abstractmethod
    def compute_survival(
        self,
        patient: Patient,
    ) -> typing.Optional[Survival]:
        pass


class SurvivalAnalysisResult:

    def __init__(
        self,
        gt_predicate: GenotypePolyPredicate,
        survival: pd.DataFrame,
        pval: float,
    ):
        assert isinstance(gt_predicate, GenotypePolyPredicate)
        self._gt_predicate = gt_predicate

        assert isinstance(survival, pd.DataFrame) and all(
            col in survival for col in ("genotype", "survival")
        )
        self._survival = survival

        if isinstance(pval, float) and math.isfinite(pval) and 0. <= pval <= 1.:
            self._pval = float(pval)
        else:
            raise ValueError(f"`p_val` must be a finite float in range [0, 1] but it was {pval}")

    @property
    def gt_predicate(self) -> GenotypePolyPredicate:
        return self._gt_predicate

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

        The DataFrame index includes the individual IDs (`patient_id`), and then there are 2 columns
        with the `genotype` group id (:attr:`~gpsea.analysis.predicate.PatientCategory.cat_id`)
        and the `survival` encoded as :class:`~gpsea.analysis.tempo.Survival` object.
        A `genotype` value may be missing (`None`) if the individual cannot be assigned
        into a genotype category.
        Similarly, a `survival` may be `None` if computing the survival is impossible for
        the individual in question.
        """
        return self._survival
    
    def complete_rows(self) -> pd.DataFrame:
        # TODO: add a convenience method for getting complete rows
        return pd.DataFrame()
    
    @property
    def pval(self) -> float:
        """
        Get the p value of the test.
        """
        return self._pval


class SurvivalAnalysis:

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
            survival=data,
            pval=pval,
        )
