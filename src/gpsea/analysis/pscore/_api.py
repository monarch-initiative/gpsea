import abc
import typing

import pandas as pd

from gpsea.model import Patient
from ..predicate.genotype import GenotypePolyPredicate
from .stats import PhenotypeScoreStatistic


class PhenotypeScorer(metaclass=abc.ABCMeta):
    """
    `PhenotypeScorer` assigns the patient with a phenotype score.
    """

    def score(self, patient: Patient) -> float:
        """
        Compute the score for the `patient`.
        """
        pass


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
        return self._genotype_phenotype_scores

    @property
    def pval(self) -> float:
        return self._pval

    def plot_boxplots(
        self,
        gt_predicate: GenotypePolyPredicate,
        ax,
        colors=["darksalmon", "honeydew"],
    ):
        # skip the patients with unassigned genotype group
        not_na_gts = self._genotype_phenotype_scores["genotype"].notna()
        data = self._genotype_phenotype_scores.loc[not_na_gts]
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
