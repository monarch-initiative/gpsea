import pytest
import ppktstore
import hpotk

from gpsea.model import Cohort, VariantEffect
from gpsea.preprocessing import configure_caching_cohort_creator, CohortCreator, load_phenopackets
from gpsea.analysis.pcats import HpoTermAnalysis
from gpsea.analysis.mtc_filter import HpoMtcFilter
from gpsea.analysis.pcats.stats import FisherExactTest
from gpsea.analysis.predicate.genotype import VariantPredicates, groups_predicate
from gpsea.analysis.predicate.phenotype import prepare_predicates_for_terms_of_interest


class TestTutorial:

    @pytest.fixture(scope='class')
    def cohort_creator(
        self,
        hpo: hpotk.MinimalOntology,
    ) -> CohortCreator:
        return configure_caching_cohort_creator(
            hpo=hpo,
        )

    @pytest.fixture(scope='class')
    def cohort(
        self,
        cohort_creator: CohortCreator,
    ) -> Cohort:
        registry = ppktstore.registry.configure_phenopacket_registry()
        with registry.open_phenopacket_store('0.1.18') as ps:
            phenopackets = tuple(ps.iter_cohort_phenopackets('TBX5'))

        cohort, _ = load_phenopackets(
            phenopackets=phenopackets,
            cohort_creator=cohort_creator,
        )
        return cohort

    @pytest.fixture
    def mtc_filter(
        self,
        hpo: hpotk.MinimalOntology,
    ) -> HpoMtcFilter:
        return HpoMtcFilter.default_filter(
            hpo=hpo,
            term_frequency_threshold=0.2,
        )

    @pytest.fixture
    def hpo_term_analysis(
        self,
        mtc_filter,
    ) -> HpoTermAnalysis:
        return HpoTermAnalysis(
            count_statistic=FisherExactTest(),
            mtc_filter=mtc_filter,
            mtc_correction='fdr_bh',
            mtc_alpha=0.05,
        )

    @pytest.mark.skip('Just for manual debugging')
    def test_compare_genotype_vs_phenotype(
        self,
        hpo: hpotk.MinimalOntology,
        cohort: Cohort,
        hpo_term_analysis: HpoTermAnalysis,
    ):
        tx_id = 'NM_181486.4'

        gt_predicate = groups_predicate(
            predicates=(
                VariantPredicates.variant_effect(VariantEffect.MISSENSE_VARIANT, tx_id),
                VariantPredicates.variant_effect(VariantEffect.FRAMESHIFT_VARIANT, tx_id),
            ),
            group_names=('Missense', 'Frameshift',),
        )
        pheno_predicates = prepare_predicates_for_terms_of_interest(
            cohort=cohort,
            hpo=hpo,
            min_n_of_patients_with_term=2,
        )

        result = hpo_term_analysis.compare_genotype_vs_phenotypes(
            cohort=cohort,
            gt_predicate=gt_predicate,
            pheno_predicates=pheno_predicates,
        )

        assert result is not None
