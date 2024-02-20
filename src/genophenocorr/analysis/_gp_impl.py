import logging
import typing

import hpotk
import pandas as pd

from statsmodels.stats import multitest

from genophenocorr.model import Patient, Cohort, VariantEffect, FeatureType
from genophenocorr.preprocessing import ProteinMetadataService
from .predicate import BooleanPredicate, PolyPredicate
from .predicate.genotype import VariantEffectPredicate, VariantPredicate, ExonPredicate, ProtFeatureTypePredicate, ProtFeaturePredicate
from .predicate.genotype import VariantEffectsPredicate, VariantsPredicate, ExonsPredicate, ProtFeaturesPredicate, ProtFeatureTypesPredicate
from .predicate.phenotype import PropagatingPhenotypeBooleanPredicateFactory, PhenotypePredicateFactory

from ._api import CohortAnalysis, GenotypePhenotypeAnalysisResult
from ._filter import PhenotypeFilter
from ._stats import run_fisher_exact, run_recessive_fisher_exact


class GpCohortAnalysis(CohortAnalysis):

    def __init__(self, cohort: Cohort,
                 hpo: hpotk.MinimalOntology,
                 protein_service: ProteinMetadataService,
                 phenotype_filter: PhenotypeFilter,
                 missing_implies_excluded: bool = False,
                 include_sv: bool = False,
                 p_val_correction: typing.Optional[str] = None,
                 ):
        if not isinstance(cohort, Cohort):
            raise ValueError(f"cohort must be type Cohort but was type {type(cohort)}")

        self._logger = logging.getLogger(__name__)
        self._hpo = hpotk.util.validate_instance(hpo, hpotk.MinimalOntology, 'hpo')
        self._phenotype_predicate_factory = PropagatingPhenotypeBooleanPredicateFactory(self._hpo,
                                                                                        missing_implies_excluded)
        self._protein_service = protein_service
        self._phenotype_filter = hpotk.util.validate_instance(phenotype_filter, PhenotypeFilter, 'phenotype_filter')
        self._correction = p_val_correction

        self._patient_list = list(cohort.all_patients) \
            if include_sv \
            else [pat for pat in cohort.all_patients if not all(var.variant_coordinates.is_structural() for var in pat.variants)]
        if len(self._patient_list) == 0:
            raise ValueError('No patients left for analysis!')

        self._hpo_terms_of_interest = self._phenotype_filter.filter_features(self._patient_list)

    def compare_by_variant_effect(self, effect: VariantEffect, tx_id: str) -> GenotypePhenotypeAnalysisResult:
        predicate = VariantEffectPredicate(tx_id, effect)
        return self._apply_boolean_predicate(predicate)

    def compare_by_variant_key(self, variant_key: str) -> GenotypePhenotypeAnalysisResult:
        predicate = VariantPredicate(variant_key)
        return self._apply_boolean_predicate(predicate)

    def compare_by_exon(self, exon_number: int, tx_id: str) -> GenotypePhenotypeAnalysisResult:
        predicate = ExonPredicate(tx_id, exon_number)
        return self._apply_boolean_predicate(predicate)

    def compare_by_protein_feature_type(self, feature_type: FeatureType, tx_id: str) -> GenotypePhenotypeAnalysisResult:
        predicate = ProtFeatureTypePredicate(tx_id, feature_type, self._protein_service)
        return self._apply_boolean_predicate(predicate)

    def compare_by_protein_feature(self, feature: str, tx_id: str) -> GenotypePhenotypeAnalysisResult:
        predicate = ProtFeaturePredicate(tx_id, feature, self._protein_service)
        return self._apply_boolean_predicate(predicate)

    def compare_by_variant_effects(self, effect1: VariantEffect, effect2: VariantEffect, tx_id: str) -> GenotypePhenotypeAnalysisResult:
        predicate = VariantEffectsPredicate(tx_id, effect1, effect2)
        return self._apply_poly_predicate(predicate)

    def compare_by_variant_keys(self, variant_key1: str, variant_key2: str) -> GenotypePhenotypeAnalysisResult:
        predicate = VariantsPredicate(variant_key1, variant_key2)
        return self._apply_poly_predicate(predicate)

    def compare_by_exons(self, exon1_number: int, exon2_number: int, tx_id: str) -> GenotypePhenotypeAnalysisResult:
        predicate = ExonsPredicate(tx_id, exon1_number, exon2_number)
        return self._apply_poly_predicate(predicate)

    def compare_by_protein_feature_types(self, feature_type1: FeatureType, feature_type2: FeatureType, tx_id: str) -> GenotypePhenotypeAnalysisResult:
        predicate = ProtFeatureTypesPredicate(tx_id, feature_type1, feature_type2, self._protein_service)
        return self._apply_poly_predicate(predicate)

    def compare_by_protein_features(self, feature1: str, feature2: str, tx_id: str) -> GenotypePhenotypeAnalysisResult:
        predicate = ProtFeaturesPredicate(tx_id, feature1, feature2, self._protein_service)
        return self._apply_poly_predicate(predicate)

    def _apply_boolean_predicate(self, predicate: BooleanPredicate) -> GenotypePhenotypeAnalysisResult:
        assert isinstance(predicate, BooleanPredicate), f'{type(predicate)} is not an instance of `BooleanPredicate`'
        return self._run_gp_analysis(self._patient_list, self._hpo_terms_of_interest,
                                     self._phenotype_predicate_factory, predicate)

    def _apply_poly_predicate(self, predicate: PolyPredicate) -> GenotypePhenotypeAnalysisResult:
        assert isinstance(predicate, PolyPredicate), f'{type(predicate)} is not an instance of `PolyPredicate`'
        return self._run_gp_analysis(self._patient_list, self._hpo_terms_of_interest,
                                     self._phenotype_predicate_factory, predicate)

    def _run_gp_analysis(self, patients: typing.Iterable[Patient],
                         phenotypes_of_interest: typing.Iterable[hpotk.TermId],
                         pheno_predicate_factory: PhenotypePredicateFactory,
                         geno_predicate: PolyPredicate) -> GenotypePhenotypeAnalysisResult:
        # 1) Count the patients
        n_usable, all_counts = self._count_patients(patients, phenotypes_of_interest,
                                                    pheno_predicate_factory, geno_predicate)

        # 2) Statistical tests
        pvals_idx = pd.Index(phenotypes_of_interest, name='p_val')
        pvals = pd.Series(float('nan'), index=pvals_idx, name='p value')
        for pf in phenotypes_of_interest:
            counts = all_counts.loc[pf]
            # TODO - this is where we must fail unless we have the contingency table of the right size!
            if counts.shape == (2, 2):
                pvals[pf] = run_fisher_exact(counts)
            elif counts.shape == (3, 2):
                pvals[pf] = run_recessive_fisher_exact(counts)
            else:
                raise ValueError(f"Invalid number of categories. A {counts.shape} table was created. Only (2, 2) and (3, 2) are valid sizes.")

        # 3) Multiple correction
        if self._correction is not None:
            _, pvals_corrected, _, _ = multitest.multipletests(pvals, alpha=.05, method=self._correction)
            corrected_idx = pd.Index(phenotypes_of_interest, name='p_val_corrected')
            corrected_pvals_series = pd.Series(data=pvals_corrected, index=corrected_idx,
                                               name='Corrected p value')
        else:
            corrected_pvals_series = None

        return GenotypePhenotypeAnalysisResult(n_usable, all_counts, pvals, corrected_pvals_series,
                                               pheno_predicate_factory.get_categories(), geno_predicate.get_question())

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
