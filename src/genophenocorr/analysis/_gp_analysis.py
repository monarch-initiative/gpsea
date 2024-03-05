import abc
import typing

import hpotk

import pandas as pd

from statsmodels.stats import multitest

from genophenocorr.model import Patient

from .predicate import PolyPredicate
from .predicate.phenotype import PhenotypePredicateFactory

from ._api import GenotypePhenotypeAnalysisResult
from ._stats import run_fisher_exact, run_recessive_fisher_exact


class GPAnalyzer(metaclass=abc.ABCMeta):
    """
    `GPAnalyzer` calculates p values for genotype-phenotype correlation of phenotypic features of interest.
    """

    def __init__(self,
                 pheno_predicate_factory: PhenotypePredicateFactory,
                 ):
        self._pheno_predicate_factory = hpotk.util.validate_instance(
            pheno_predicate_factory, PhenotypePredicateFactory, 'pheno_predicate_factory')

    @abc.abstractmethod
    def analyze(self,
                patients: typing.Iterable[Patient],
                phenotypic_features: typing.Iterable[hpotk.TermId],
                predicate: PolyPredicate,
                ) -> GenotypePhenotypeAnalysisResult:
        """
        Test for association between `phenotypic_features` of interest groups or `patients` determined
        by the `predicate`.

        :return: :class:`GenotypePhenotypeAnalysisResult` object with the results.
        """
        pass

    @staticmethod
    def _count_patients(patients: typing.Iterable[Patient],
                        phenotypes_of_interest: typing.Iterable[hpotk.TermId],
                        pheno_predicate_factory: PhenotypePredicateFactory,
                        geno_predicate: PolyPredicate) -> typing.Tuple[pd.Series, pd.DataFrame]:
        row_idx = pd.MultiIndex.from_product(
            iterables=(phenotypes_of_interest, pheno_predicate_factory.get_categories()),
            names=('phenotypic_feature', 'category')
        )
        col_idx = pd.Index(geno_predicate.get_categories(), name=geno_predicate.get_question())
        counts = pd.DataFrame(data=0, index=row_idx, columns=col_idx)
        n_usable_patients = pd.Series(data=0, index=pd.Index(phenotypes_of_interest))

        # Apply genotype and phenotype predicates
        for pf in phenotypes_of_interest:
            pheno_predicate = pheno_predicate_factory.get_predicate(pf)

            for patient in patients:
                pheno_cat = pheno_predicate.test(patient)
                geno_cat = geno_predicate.test(patient)

                if pheno_cat is not None and geno_cat is not None:
                    counts.loc[(pf, pheno_cat), geno_cat] += 1
                    n_usable_patients[pf] += 1

        return n_usable_patients, counts


class FisherExactAnalyzer(GPAnalyzer):
    """
    `FisherExactAnalyzer` uses Fisher's exact test to calculate p value for phenotypic features of interest.

    Following the test, the code applies one of the multiple testing corrections provided by
    :func:`statsmodels.stats.multitest.multipletests` function at given `mtc_alpha`.

    :param p_val_correction: a `str` with name of the multiple testing correction method or `None` if no correction should be applied.
    :param mtc_alpha: a `float` in range :math:`(0, 1]` with the multiple testing correction alpha value.
    """

    def __init__(self,
                 pheno_predicate_factory: PhenotypePredicateFactory,
                 p_val_correction: typing.Optional[str] = None,
                 mtc_alpha: float = .05,
                 ):
        super().__init__(pheno_predicate_factory)
        self._correction = p_val_correction
        self._mtc_alpha = mtc_alpha

    def analyze(
            self,
            patients: typing.Iterable[Patient],
            phenotypic_features: typing.Iterable[hpotk.TermId],
            predicate: PolyPredicate,
    ) -> GenotypePhenotypeAnalysisResult:
        # 1) Count the patients
        n_usable, all_counts = GPAnalyzer._count_patients(
            patients, phenotypic_features,
            self._pheno_predicate_factory, predicate,
        )

        # 2) Statistical tests
        pvals_idx = pd.Index(phenotypic_features, name='p_val')
        pvals = pd.Series(float('nan'), index=pvals_idx, name='p value')
        for pf in phenotypic_features:
            counts = all_counts.loc[pf]
            # TODO - this is where we must fail unless we have the contingency table of the right size!
            if counts.shape == (2, 2):
                pvals[pf] = run_fisher_exact(counts)
            elif counts.shape == (3, 2):
                pvals[pf] = run_recessive_fisher_exact(counts)
            else:
                raise ValueError(
                    "Invalid number of categories. "
                    f"A {counts.shape} table was created. Only (2, 2) and (3, 2) are valid sizes."
                )

        # 3) Multiple correction
        if self._correction is not None:
            _, pvals_corrected, _, _ = multitest.multipletests(pvals, alpha=self._mtc_alpha, method=self._correction)
            corrected_idx = pd.Index(phenotypic_features, name='p_val_corrected')
            corrected_pvals_series = pd.Series(data=pvals_corrected, index=corrected_idx,
                                               name='Corrected p value')
        else:
            corrected_pvals_series = None

        # 4) Wrap up
        return GenotypePhenotypeAnalysisResult(
            n_usable=n_usable,
            all_counts=all_counts,
            pvals=pvals,
            corrected_pvals=corrected_pvals_series,
            phenotype_categories=self._pheno_predicate_factory.get_categories(),
            question=predicate.get_question(),
        )
