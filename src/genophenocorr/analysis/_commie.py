import logging
import typing
from collections import Counter

import hpotk
import pandas as pd

from statsmodels.stats import multitest

from genophenocorr.model import Patient, Cohort, VariantEffect, FeatureType
from .predicate import BooleanPredicate, PolyPredicate
from .predicate.genotype import VariantEffectPredicate, VariantPredicate, ExonPredicate, ProtFeatureTypePredicate, ProtFeaturePredicate, VariantsPredicate
from .predicate.phenotype import PropagatingPhenotypeBooleanPredicateFactory, PhenotypePredicateFactory

from ._api import CohortAnalysis, GenotypePhenotypeAnalysisResult
from ._stats import run_fisher_exact


def _filter_rare_phenotypes_using_hierarchy(patients: typing.Collection[Patient],
                                            min_perc_patients_w_hpo: float,
                                            hpo: hpotk.GraphAware,
                                            missing_implies_excluded: bool = False) -> typing.Sequence[hpotk.TermId]:
    """
    We want to get a corpus of HPO terms that are present in at least certain fraction of the subjects.
    We will get the IDs of the present features only.

    :param min_perc_patients_w_hpo: the minimum fraction of patients that *must* have the feature.
      The `float` must be in range :math:`(0,1]`.
    """
    # TODO - review the logic with Lauren. Should we also do annotation propagation here?
    #  What is the meaning of _include_unmeasured in this context?
    present_count = Counter()
    excluded_count = Counter()

    for patient in patients:
        for pf in patient.phenotypes:
            if pf.is_observed:
                # A present phenotypic feature must be counted in.
                present_count[pf.identifier] += 1
                # implies presence of its ancestors.
                for anc in hpo.graph.get_ancestors(pf):
                    present_count[anc] += 1
            else:
                # An excluded phenotypic feature
                excluded_count[pf.identifier] += 1
                for desc in hpo.graph.get_descendants(pf):
                    # implies exclusion of its descendants.
                    excluded_count[desc] += 1

    if missing_implies_excluded:
        # We must do another pass to ensure that we account the missing features as excluded.
        # `present_count` has all phenotypic features that were observed as present among the cohort members.
        for pf in present_count:
            for patient in patients:
                pf_can_be_inferred_as_excluded_for_patient = False
                for phenotype in patient.phenotypes:
                    if pf == phenotype.identifier:
                        # Patient is explicitly annotated with `pf`. No inference for this patient!
                        pf_can_be_inferred_as_excluded_for_patient = False
                        break

                    if phenotype.is_observed and hpo.graph.is_ancestor_of(pf, phenotype):
                        # The `pf` is implicitly observed in the patient. No inference for this patient!
                        pf_can_be_inferred_as_excluded_for_patient = False
                        break
                    elif not phenotype.is_observed and hpo.graph.is_descendant_of(pf, phenotype):
                        # The `pf` is implicitly excluded in the patient. No inference for this patient!
                        pf_can_be_inferred_as_excluded_for_patient = False
                        break
                    else:
                        # The `pf` is observed or excluded neither implicitly nor explicitly.
                        # We can infer that it is excluded!
                        pf_can_be_inferred_as_excluded_for_patient = True

                if pf_can_be_inferred_as_excluded_for_patient:
                    excluded_count[pf] += 1
                    for desc in hpo.graph.get_descendants(pf):
                        excluded_count[desc] += 1

    total_count = Counter()
    for term_id, count in present_count.items():
        total_count[term_id] += count

    for term_id, count in excluded_count.items():
        total_count[term_id] += count

    final_hpo = []
    for term_id, n_present in present_count.items():
        n_all = total_count[term_id]
        ratio = n_present / n_all
        if ratio >= min_perc_patients_w_hpo:
            final_hpo.append(term_id)

    if len(final_hpo) == 0:
        raise ValueError(f"No HPO terms found in over {min_perc_patients_w_hpo * 100}% of patients.")

    return final_hpo



class CommunistCohortAnalysis(CohortAnalysis):

    def __init__(self, cohort: Cohort,
                 hpo: hpotk.MinimalOntology,
                 missing_implies_excluded: bool = False,
                 include_sv: bool = False,
                 p_val_correction: typing.Optional[str] = None,
                 min_perc_patients_w_hpo: typing.Union[float, int] = .1):
        if not isinstance(cohort, Cohort):
            raise ValueError(f"cohort must be type Cohort but was type {type(cohort)}")

        self._logger = logging.getLogger(__name__)
        self._hpo = hpotk.util.validate_instance(hpo, hpotk.MinimalOntology, 'hpo')
        self._phenotype_predicate_factory = PropagatingPhenotypeBooleanPredicateFactory(self._hpo,
                                                                                        missing_implies_excluded)
        self._correction = p_val_correction
        self._patient_list = list(cohort.all_patients) \
            if include_sv \
            else [pat for pat in cohort.all_patients if not all(var.variant_coordinates.is_structural() for var in pat.variants)]
        if len(self._patient_list) == 0:
            raise ValueError('No patients left for analysis!')
        self._percent_pats = self._check_min_perc_patients_w_hpo(min_perc_patients_w_hpo, len(self._patient_list))

        # As of now, `self._testing_hpo_terms` is a Sequence[TermId] with terms that were *present* in >=1 cohort member.
        self._testing_hpo_terms = _filter_rare_phenotypes_using_hierarchy(self._patient_list, self._percent_pats,
                                                                          self._hpo, missing_implies_excluded)

        self._patients_by_hpo = self._group_patients_by_hpo(self._testing_hpo_terms, self._patient_list,
                                                            self._hpo, missing_implies_excluded)

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
        predicate = ProtFeatureTypePredicate(transcript_id=tx_id, feature_type=feature_type)
        return self._apply_boolean_predicate(predicate)

    def compare_by_protein_feature(self, feature: str, tx_id: str) -> GenotypePhenotypeAnalysisResult:
        predicate = ProtFeaturePredicate(tx_id, feature)
        return self._apply_boolean_predicate(predicate)

    def compare_by_variant_keys(self, a: str, b: str) -> GenotypePhenotypeAnalysisResult:
        predicate = VariantsPredicate(a, b)
        return self._apply_poly_predicate(predicate)

    def _apply_boolean_predicate(self, predicate: BooleanPredicate) -> GenotypePhenotypeAnalysisResult:
        return self._run_gp_analysis(self._patient_list, self._testing_hpo_terms,
                                     self._phenotype_predicate_factory, predicate)

    def _apply_poly_predicate(self, predicate: PolyPredicate) -> GenotypePhenotypeAnalysisResult:
        return self._run_gp_analysis(self._patient_list, self._testing_hpo_terms,
                                     self._phenotype_predicate_factory, predicate)

    def _run_gp_analysis(self, patients: typing.Iterable[Patient],
                         phenotypic_features: typing.Iterable[hpotk.TermId],
                         pheno_predicate_factory: PhenotypePredicateFactory,
                         geno_predicate: PolyPredicate) -> GenotypePhenotypeAnalysisResult:
        # 1) Count the patients
        n_usable, all_counts = self._count_patients(patients, phenotypic_features,
                                                    pheno_predicate_factory, geno_predicate)

        # 2) Statistical tests
        pvals_idx = pd.Index(phenotypic_features, name='p_val')
        pvals = pd.Series(float('nan'), index=pvals_idx, name='p value')
        for pf in phenotypic_features:
            counts = all_counts.loc[pf]
            # TODO - this is where we must fail unless we have the contingency table of the right size!
            pvals[pf] = run_fisher_exact(counts)

        # 3) Multiple correction
        if self._correction is not None:
            _, pvals_corrected, _, _ = multitest.multipletests(pvals, alpha=.05, method=self._correction)
            corrected_idx = pd.Index(phenotypic_features, name='p_val_corrected')
            corrected_pvals_series = pd.Series(data=pvals_corrected, index=corrected_idx,
                                               name='Corrected p value')
        else:
            corrected_pvals_series = None

        return GenotypePhenotypeAnalysisResult(n_usable, all_counts, pvals, corrected_pvals_series,
                                               pheno_predicate_factory.get_categories(), geno_predicate.get_question())

    @staticmethod
    def _count_patients(patients: typing.Iterable[Patient],
                        phenotypic_features: typing.Iterable[hpotk.TermId],
                        pheno_predicate_factory: PhenotypePredicateFactory,
                        geno_predicate: PolyPredicate) -> typing.Tuple[pd.Series, pd.DataFrame]:
        row_idx = pd.MultiIndex.from_product(
            iterables=(phenotypic_features, pheno_predicate_factory.get_categories()),
            names=('phenotypic_feature', 'category')
        )
        col_idx = pd.Index(geno_predicate.get_categories(), name=geno_predicate.get_question())
        counts = pd.DataFrame(data=0, index=row_idx, columns=col_idx)
        n_usable_patients = pd.Series(data=0, index=pd.Index(phenotypic_features))

        # Apply genotype and phenotype predicates
        for pf in phenotypic_features:
            pheno_predicate = pheno_predicate_factory.get_predicate(pf)

            for patient in patients:
                pheno_cat = pheno_predicate.test(patient)
                geno_cat = geno_predicate.test(patient)

                if pheno_cat is not None and geno_cat is not None:
                    counts.loc[(pf, pheno_cat), geno_cat] += 1
                    n_usable_patients[pf] += 1

        return n_usable_patients, counts
