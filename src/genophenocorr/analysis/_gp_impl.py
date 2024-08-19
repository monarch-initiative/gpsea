import logging
import typing

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
from .predicate.phenotype import PhenotypePolyPredicate, P, PropagatingPhenotypePredicate, DiseasePresencePredicate, CountingPhenotypeScorer

from ._api import CohortAnalysis, GenotypePhenotypeAnalysisResult, PhenotypeScoreAnalysisResult
from ._filter import PhenotypeFilter
from ._gp_analysis import GPAnalyzer


class GpCohortAnalysis(CohortAnalysis):

    def __init__(
        self, cohort: Cohort,
        hpo: hpotk.MinimalOntology,
        protein_service: ProteinMetadataService,
        phenotype_filter: PhenotypeFilter,
        gp_analyzer: GPAnalyzer,
        missing_implies_excluded: bool,
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
        self._phenotype_filter = hpotk.util.validate_instance(phenotype_filter, PhenotypeFilter, 'phenotype_filter')
        self._gp_analyzer = hpotk.util.validate_instance(gp_analyzer, GPAnalyzer, 'gp_analyzer')

        self._patient_list = list(cohort.all_patients) \
            if include_sv \
            else [pat for pat in cohort.all_patients if not all(var.variant_info.is_structural() for var in pat.variants)]
        if len(self._patient_list) == 0:
            raise ValueError('No patients left for analysis!')

        self._hpo_terms_of_interest = self._phenotype_filter.filter_features(self._patient_list)
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
        return self._apply_poly_predicate_on_hpo_terms(bool_predicate)

    def compare_hpo_vs_recessive_genotype(
        self,
        predicate: VariantPredicate,
    ) -> GenotypePhenotypeAnalysisResult:
        """
        Bin patients according to a presence of zero, one, or two alleles that matche the `predicate`
        and test for genotype-phenotype correlations.
        """
        rec_predicate = wrap_as_recessive_predicate(predicate)
        return self._apply_poly_predicate_on_hpo_terms(rec_predicate)

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
        return self._apply_poly_predicate_on_hpo_terms(predicate)

    def compare_disease_vs_genotype(
        self,
        predicate: VariantPredicate,
        disease_ids: typing.Optional[typing.Sequence[typing.Union[str, hpotk.TermId]]] = None,
    ) -> GenotypePhenotypeAnalysisResult:
        pheno_predicates = self._prepare_disease_predicates(disease_ids)

        # This can be updated to any genotype poly predicate in future, if necessary.
        genotype_predicate = wrap_as_boolean_predicate(predicate)
        return self._apply_poly_predicate(pheno_predicates, genotype_predicate)

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
            data.loc[patient.patient_id, 'genotype'] = None if gt_cat is None else gt_cat.category.name
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

    def _apply_poly_predicate_on_hpo_terms(
            self,
            gt_predicate: GenotypePolyPredicate,
    ) -> GenotypePhenotypeAnalysisResult:
        assert isinstance(gt_predicate, GenotypePolyPredicate), \
            f'{type(gt_predicate)} is not an instance of `GenotypePolyPredicate`'

        pheno_predicates = self._prepare_phenotype_predicates()
        return self._apply_poly_predicate(pheno_predicates, gt_predicate)

    # Make public, eventually convenience functions written
    def _apply_poly_predicate(
            self,
            pheno_predicates: typing.Iterable[PhenotypePolyPredicate[P]],
            gt_predicate: GenotypePolyPredicate,
    ):
        return self._gp_analyzer.analyze(
            patients=self._patient_list,
            pheno_predicates=pheno_predicates,
            gt_predicate=gt_predicate,
        )

    def _prepare_phenotype_predicates(self) -> typing.Sequence[PhenotypePolyPredicate[P]]:
        return tuple(
            PropagatingPhenotypePredicate(
                hpo=self._hpo,
                query=query,
            )
            for query in self._hpo_terms_of_interest)
