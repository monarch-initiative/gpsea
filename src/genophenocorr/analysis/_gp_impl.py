import logging
import typing

from collections import Counter

import hpotk
import pandas as pd

from scipy.stats import mannwhitneyu

from genophenocorr.model import Cohort, Patient
from genophenocorr.preprocessing import ProteinMetadataService
from .predicate import GenotypePolyPredicate
from .predicate.genotype import VariantPredicate
from .predicate.genotype import boolean_predicate as wrap_as_boolean_predicate
from .predicate.genotype import groups_predicate as wrap_as_groups_predicate
from .predicate.genotype import recessive_predicate as wrap_as_recessive_predicate
from .predicate.phenotype import PhenotypePolyPredicate, P, PropagatingPhenotypePredicate, DiseasePresencePredicate

from ._api import CohortAnalysis, GenotypePhenotypeAnalysisResult, PhenotypeScoreAnalysisResult
from ._gp_analysis import GPAnalyzer


class GpCohortAnalysis(CohortAnalysis):

    def __init__(
        self, cohort: Cohort,
        hpo: hpotk.MinimalOntology,
        protein_service: ProteinMetadataService,
        gp_analyzer: GPAnalyzer,
        missing_implies_excluded: bool,
        min_n_of_patients_with_term: int,
        include_sv: bool = False,
    ):
        super().__init__(
            hpo,
            protein_service,
        )
        if not isinstance(cohort, Cohort):
            raise ValueError(f"cohort must be type Cohort but was type {type(cohort)}")

        self._logger = logging.getLogger(__name__)
        self._cohort = cohort
        # self._phenotype_filter = hpotk.util.validate_instance(phenotype_filter, PhenotypeFilter, 'phenotype_filter')
        self._gp_analyzer = hpotk.util.validate_instance(gp_analyzer, GPAnalyzer, 'gp_analyzer')

        self._patient_list = list(cohort.all_patients) \
            if include_sv \
            else [pat for pat in cohort.all_patients if not all(var.variant_info.is_structural() for var in pat.variants)]
        if len(self._patient_list) == 0:
            raise ValueError('No patients left for analysis!')

        self._hpo_terms_of_interest = prepare_hpo_terms_of_interest(
            hpo=self._hpo,
            patients=self._patient_list,
            missing_implies_excluded=missing_implies_excluded,
            min_n_of_patients_with_term=min_n_of_patients_with_term,
        )
        self._missing_implies_excluded = missing_implies_excluded

    def compare_hpo_vs_genotype(
            self,
            predicate: VariantPredicate,
    ) -> GenotypePhenotypeAnalysisResult:
        """
        Bin patients according to a presence of at least one allele that matches `predicate`
        and test for genotype-phenotype correlations.
        """
        assert isinstance(predicate, VariantPredicate), \
            f'{type(predicate)} is not an instance of `VariantPredicate`'
        bool_predicate = wrap_as_boolean_predicate(predicate) 
        return self.compare_genotype_vs_cohort_phenotypes(bool_predicate)

    def compare_hpo_vs_recessive_genotype(
        self,
        predicate: VariantPredicate,
    ) -> GenotypePhenotypeAnalysisResult:
        """
        Bin patients according to a presence of zero, one, or two alleles that matche the `predicate`
        and test for genotype-phenotype correlations.
        """
        rec_predicate = wrap_as_recessive_predicate(predicate)
        return self.compare_genotype_vs_cohort_phenotypes(rec_predicate)

    def compare_hpo_vs_genotype_groups(
        self,
        predicates: typing.Iterable[VariantPredicate],
        group_names: typing.Iterable[str],
    ) -> GenotypePhenotypeAnalysisResult:
        """
        Bin patients according to a presence of at least one allele that matches
        any of the provided `predicates` and test for genotype-phenotype correlations
        between the groups.

        Note, the patients that pass testing by >1 genotype predicate are *OMITTED* from the analysis!
        """
        predicate = wrap_as_groups_predicate(
            predicates=predicates,
            group_names=group_names,
        )
        return self.compare_genotype_vs_cohort_phenotypes(predicate)

    def compare_disease_vs_genotype(
        self,
        predicate: VariantPredicate,
        disease_ids: typing.Optional[typing.Sequence[typing.Union[str, hpotk.TermId]]] = None,
    ) -> GenotypePhenotypeAnalysisResult:
        pheno_predicates = self._prepare_disease_predicates(disease_ids)

        # This can be updated to any genotype poly predicate in future, if necessary.
        genotype_predicate = wrap_as_boolean_predicate(predicate)
        return self.compare_genotype_vs_phenotypes(pheno_predicates, genotype_predicate)

    def _prepare_disease_predicates(
        self,
        disease_ids: typing.Optional[typing.Sequence[typing.Union[str, hpotk.TermId]]],
    ) -> typing.Sequence[DiseasePresencePredicate]:
        testing_diseases = []
        if disease_ids is None:
            disease_ids = [dis.identifier for dis in self._cohort.all_diseases()]
        if len(disease_ids) < 1:
            raise ValueError("No diseases available for testing.")
        for disease_id in disease_ids:
            if isinstance(disease_id, str):
                testing_diseases.append(hpotk.TermId.from_curie(disease_id))
            elif isinstance(disease_id, hpotk.TermId):
                testing_diseases.append(disease_id)
            else:
                raise ValueError(f'{disease_id} must be a `str` or a `hpotk.TermId`')
        pheno_predicates = []
        for disease in testing_diseases:
            pheno_predicates.append(DiseasePresencePredicate(disease))
        return pheno_predicates

    def compare_genotype_vs_phenotype_score(
        self,
        gt_predicate: GenotypePolyPredicate,
        phenotype_scorer: typing.Callable[[Patient,], float],
    ) -> PhenotypeScoreAnalysisResult:
        idx = pd.Index(tuple(p.patient_id for p in self._patient_list), name='patient_id')
        data = pd.DataFrame(
            None,
            index=idx,
            columns=['genotype', 'phenotype'],
        )

        # Apply the predicates on the patients
        for patient in self._patient_list:
            gt_cat = gt_predicate.test(patient)
            data.loc[patient.patient_id, 'genotype'] = None if gt_cat is None else gt_cat.category.cat_id
            data.loc[patient.patient_id, 'phenotype'] = phenotype_scorer(patient)

        # To improve the determinism
        data.sort_index(inplace=True)

        # Sort by PatientCategory.cat_id and unpack.
        # For now, we only allow to have up to 2 groups.
        x_key, y_key = sorted(data['genotype'].dropna().unique())
        x = data.loc[data['genotype'] == x_key, 'phenotype'].to_numpy(dtype=float)
        y = data.loc[data['genotype'] == y_key, 'phenotype'].to_numpy(dtype=float)
        result = mannwhitneyu(
            x=x,
            y=y,
            alternative='two-sided',
        )

        return PhenotypeScoreAnalysisResult(
            genotype_phenotype_scores=data,
            p_value=float(result.pvalue),
        )

    def compare_genotype_vs_cohort_phenotypes(
        self,
        gt_predicate: GenotypePolyPredicate,
    ) -> GenotypePhenotypeAnalysisResult:
        assert isinstance(gt_predicate, GenotypePolyPredicate), \
            f'{type(gt_predicate)} is not an instance of `GenotypePolyPredicate`'

        pheno_predicates = self._prepare_phenotype_predicates()
        return self.compare_genotype_vs_phenotypes(gt_predicate, pheno_predicates)

    def _prepare_phenotype_predicates(self) -> typing.Sequence[PhenotypePolyPredicate[P]]:
        return tuple(
            PropagatingPhenotypePredicate(
                hpo=self._hpo,
                query=query,
            )
            for query in self._hpo_terms_of_interest)

    def compare_genotype_vs_phenotypes(
        self,
        gt_predicate: GenotypePolyPredicate,
        pheno_predicates: typing.Iterable[PhenotypePolyPredicate[P]],
    ) -> GenotypePhenotypeAnalysisResult:
        return self._gp_analyzer.analyze(
            patients=self._patient_list,
            pheno_predicates=pheno_predicates,
            gt_predicate=gt_predicate,
        )


def prepare_hpo_terms_of_interest(
    hpo: hpotk.MinimalOntology,
    patients: typing.Iterable[Patient],
    missing_implies_excluded: bool,
    min_n_of_patients_with_term: int,
) -> typing.Iterable[hpotk.TermId]:
    present_count = Counter()
    excluded_count = Counter()

    for patient in patients:
        for pf in patient.phenotypes:
            if pf.is_present:
                # A present phenotypic feature must be counted in.
                present_count[pf.identifier] += 1
                # and it also implies presence of its ancestors.
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

                    if phenotype.is_present and hpo.graph.is_ancestor_of(pf, phenotype):
                        # The `pf` is implicitly observed in the patient. No inference for this patient!
                        pf_can_be_inferred_as_excluded_for_patient = False
                        break
                    elif phenotype.is_excluded and hpo.graph.is_descendant_of(pf, phenotype):
                        # The `pf` is implicitly excluded in the patient. No inference for this patient!
                        pf_can_be_inferred_as_excluded_for_patient = False
                        break
                    else:
                        # There is no implicit or explicit observation of `pf` for this patient.
                        # Therefore, we infer that it is excluded!
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
    for term_id in present_count:
        # Keep the term if it is mentioned at least *n* times (incl. being excluded)
        # in the cohort
        n_all = total_count[term_id]
        if n_all >= min_n_of_patients_with_term:
            final_hpo.append(term_id)

    return tuple(final_hpo)
