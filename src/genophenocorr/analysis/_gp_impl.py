import logging
import typing

import hpotk

from genophenocorr.model import Cohort, VariantEffect, FeatureType
from genophenocorr.model.genome import Region
from genophenocorr.preprocessing import ProteinMetadataService
from .predicate import GenotypePolyPredicate
from .predicate.genotype import RecessiveVariantPredicate, RecessiveProtFeaturePredicate, RecessiveExonPredicate, RecessiveProtFeatureTypePredicate, RecessiveVariantEffectPredicate, RecessiveProtRegionPredicate
from .predicate.genotype import VariantPredicate
from .predicate.genotype import boolean_predicate as wrap_as_boolean_predicate, grouping_predicate as wrap_as_grouping_predicate
from .predicate.phenotype import PhenotypePolyPredicate, P, PropagatingPhenotypePredicate, DiseasePresencePredicate

from ._api import CohortAnalysis, GenotypePhenotypeAnalysisResult
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
        super().__init__(protein_service)
        if not isinstance(cohort, Cohort):
            raise ValueError(f"cohort must be type Cohort but was type {type(cohort)}")

        self._logger = logging.getLogger(__name__)
        self._cohort = cohort
        self._hpo = hpotk.util.validate_instance(hpo, hpotk.MinimalOntology, 'hpo')
        self._phenotype_filter = hpotk.util.validate_instance(phenotype_filter, PhenotypeFilter, 'phenotype_filter')
        self._gp_analyzer = hpotk.util.validate_instance(gp_analyzer, GPAnalyzer, 'gp_analyzer')

        self._patient_list = list(cohort.all_patients) \
            if include_sv \
            else [pat for pat in cohort.all_patients if not all(var.variant_coordinates.is_structural() for var in pat.variants)]
        if len(self._patient_list) == 0:
            raise ValueError('No patients left for analysis!')

        self._hpo_terms_of_interest = self._phenotype_filter.filter_features(self._patient_list)
        self._missing_implies_excluded = missing_implies_excluded
    
    # TODO: remove
    def compare_by_recessive_variant_effect(self, effect: VariantEffect, tx_id: str) -> GenotypePhenotypeAnalysisResult:
        predicate = RecessiveVariantEffectPredicate(tx_id, effect)
        return self._apply_poly_predicate_on_hpo_terms(predicate)

    def compare_by_recessive_variant_key(self, variant_key: str) -> GenotypePhenotypeAnalysisResult:
        predicate = RecessiveVariantPredicate(variant_key)
        return self._apply_poly_predicate_on_hpo_terms(predicate)

    def compare_by_recessive_exon(self, exon_number: int, tx_id: str) -> GenotypePhenotypeAnalysisResult:
        predicate = RecessiveExonPredicate(tx_id, exon_number)
        return self._apply_poly_predicate_on_hpo_terms(predicate)

    def compare_by_recessive_protein_feature_type(self, feature_type: FeatureType, tx_id: str) -> GenotypePhenotypeAnalysisResult:
        predicate = RecessiveProtFeatureTypePredicate(tx_id, feature_type, self._protein_service)
        return self._apply_poly_predicate_on_hpo_terms(predicate)

    def compare_by_recessive_protein_feature(self, feature: str, tx_id: str) -> GenotypePhenotypeAnalysisResult:
        predicate = RecessiveProtFeaturePredicate(tx_id, feature, self._protein_service)
        return self._apply_poly_predicate_on_hpo_terms(predicate)
    
    def compare_by_recessive_protein_region(self, region:Region, tx_id:str) -> GenotypePhenotypeAnalysisResult:
        predicate = RecessiveProtRegionPredicate(tx_id, region)
        return self._apply_poly_predicate_on_hpo_terms(predicate)

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
        predicate = wrap_as_boolean_predicate(predicate) 
        return self._apply_poly_predicate_on_hpo_terms(predicate)
    
    def compare_hpo_vs_genotype_groups(
            self,
            first: VariantPredicate,
            second: VariantPredicate,
    ) -> GenotypePhenotypeAnalysisResult:
        """
        Bin patients according to a presence of at least one allele that matches `first` or `second` predicate 
        and test for genotype-phenotype correlations between the groups.

        Note, the patients that pass testing by both predicates are *OMITTED* from the analysis!
        """
        predicate = wrap_as_grouping_predicate(
            first=first,
            second=second,
        )
        return self._apply_poly_predicate_on_hpo_terms(predicate)

    def compare_disease_vs_genotype(
        self, 
        genotype_predicate: GenotypePolyPredicate, 
        disease_ids: typing.Optional[typing.Sequence[typing.Union[str, hpotk.TermId]]] = None,
    ) -> GenotypePhenotypeAnalysisResult:
        pheno_predicates = []
        testing_diseases = []
        if disease_ids is None:
            disease_ids = [dis.identifier for dis in self._cohort.all_diseases()]
        if len(disease_ids) < 1:
            raise ValueError("No diseases available for testing.")
        for dis in disease_ids:
            if type(dis) == str:
                testing_diseases.append(hpotk.TermId.from_curie(dis))
            else:
                testing_diseases.append(hpotk.util.validate_instance(dis, hpotk.TermId, 'dis'))

        for disease in testing_diseases:
            pheno_predicates.append(DiseasePresencePredicate(disease))

        return self._apply_poly_predicate(pheno_predicates, genotype_predicate)

    def _apply_poly_predicate_on_hpo_terms(
            self,
            gt_predicate: GenotypePolyPredicate,
    ) -> GenotypePhenotypeAnalysisResult:
        assert isinstance(gt_predicate, GenotypePolyPredicate), \
            f'{type(gt_predicate)} is not an instance of `GenotypePolyPredicate`'

        pheno_predicates = self._prepare_phenotype_predicates()
        return self._apply_poly_predicate(pheno_predicates, gt_predicate)

    #Make public, eventually convenience functions written
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
