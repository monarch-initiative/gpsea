import typing

import hpotk
import pytest

import pandas as pd

from genophenocorr.analysis import configure_cohort_analysis, GenotypePhenotypeAnalysisResult
from genophenocorr.analysis.predicate import PatientCategories
from genophenocorr.analysis.predicate.genotype import VariantPredicates
from genophenocorr.model import Cohort, VariantEffect


class TestCohortAnalysis:

    def test_compare_by_variant_effect(
            self,
            suox_cohort: Cohort,
            hpo: hpotk.MinimalOntology,
    ):
        pd.set_option('expand_frame_repr', False)
        cohort_analysis = configure_cohort_analysis(suox_cohort, hpo)
        results = cohort_analysis.compare_by_variant_effect(VariantEffect.MISSENSE_VARIANT, 'NM_001032386.2')
        print(results)
        summary = results.summarize(hpo, PatientCategories.YES)
        print(summary)

    def test_get_count(
            self,
            suox_cohort: Cohort,
            hpo: hpotk.MinimalOntology,
    ):
        """
        This test shows how to manipulate the results object to get the counts we need
        Let's use Arachnodactyly [HP:0001166]  Yes (Genotype NO 1/10); No (Genotype YES 13/16) as an example
        that is
        Arachnodactyly (+) MISSENSE (-) = 1
        Arachnodactyly (-) MISSENSE (-) = 9
        Arachnodactyly (+) MISSENSE (+) = 13
        Arachnodactyly (-) MISSENSE (+) = 3
        In the Toy file, we have
        Arachnodactyly TRUE, MISSENSE (snv) TRUE: A,B,D,E;G;J;M;P;Q;R;T;V;Y = 13
        Arachnodactyly FALSE, MISSENSE (snv) TRUE: C,K,N=3
        Arachnodactyly TRUE, MISSENSE (snv) FALSE: H = 1
        Arachnodactyly FALSE, MISSENSE (snv) FALSE: F,I,L,O,S;U,W,X,Z=9
        This is the output
        MISSENSE_VARIANT on NM_1234.5                         No            Yes
                                                   Count Percent  Count Percent   p value Corrected p value
        Arachnodactyly [HP:0001166]                         1/10   10.0%  13/16  81.25%  0.000781          0.020299
        """
        pd.set_option('expand_frame_repr', False)
        cohort_analysis = configure_cohort_analysis(suox_cohort, hpo)
        results = cohort_analysis.compare_by_variant_effect(VariantEffect.MISSENSE_VARIANT, 'NM_001032386.2')

        # Let's make sure we know what class we have
        assert isinstance(results, GenotypePhenotypeAnalysisResult)

        # We expect that we have YES and NO as phenotype categories
        phenotype_categories_tuple = results.phenotype_categories
        assert len(phenotype_categories_tuple) == 2
        assert PatientCategories.YES in phenotype_categories_tuple
        assert PatientCategories.NO in phenotype_categories_tuple

        # The all_counts field has the counts of the phenotypes
        all_counts = results.all_counts
        assert isinstance(all_counts, typing.Mapping)

        # We tested 74 HPO terms
        assert len(all_counts) == 74

        # The index of all_counts is a Tuple with (HPO TermId, BooleanPredicate
        # Let's test Seizure - we should have one row for each Patient Predicate
        counts = all_counts[hpotk.TermId.from_curie("HP:0001250")]

        # The YES row is Seizure YES -- according to the above, we have 11 (MISSENSE NO) and 17 (MISSENSE YES)
        assert counts.loc[PatientCategories.YES, PatientCategories.NO] == 11
        assert counts.loc[PatientCategories.YES, PatientCategories.YES] == 17

        # The NO row is Seizure NO -- according to the above, we have 0 (MISSENSE NO) and 7 (MISSENSE YES)
        assert counts.loc[PatientCategories.NO, PatientCategories.NO] == 0
        assert counts.loc[PatientCategories.NO, PatientCategories.YES] == 7

        # In total, 35 patients were categorized
        assert counts.sum().sum() == 35

    def test_get_positive_count(
            self,
            suox_cohort: Cohort,
            hpo: hpotk.MinimalOntology,
    ):
        """
        This test shows how to get the counts for positive HPO terms
        Let's use Seizure [HP:0001250]  Yes (Genotype NO 11/11); No (Genotype YES 17/24) as an example
        that is
        Seizure (+) MISSENSE (-) = 11
        Seizure (-) MISSENSE (-) = 0
        Seizure (+) MISSENSE (+) = 17
        Seizure (-) MISSENSE (+) = 7
        This means we expect 11+17=28
        See the previous test for further information
        """
        pd.set_option('expand_frame_repr', False)
        cohort_analysis = configure_cohort_analysis(suox_cohort, hpo)
        results = cohort_analysis.compare_by_variant_effect(VariantEffect.MISSENSE_VARIANT, 'NM_001032386.2')
        all_counts = results.all_counts
        # The index of all_counts is a Tuple with (HPO TermId, BooleanPredicate
        # Let's test Arachnodactyly - we should have one row for each Patient Predicate
        counts = all_counts[hpotk.TermId.from_curie("HP:0001250")]
        #total_observed_HPO = HeuristicSamplerMtcFilter.get_number_of_positive_observations(arachnodactyly_counts)
        total_observed_HPO = counts.loc[PatientCategories.YES, PatientCategories.NO] + counts.loc[PatientCategories.YES, PatientCategories.YES]
        assert 28 == total_observed_HPO

    @pytest.mark.skip('For debugging only')
    def test_compare_symptom_count_vs_genotype(
        self,
        suox_cohort: Cohort,
        hpo: hpotk.MinimalOntology,
    ):
        cohort_analysis = configure_cohort_analysis(suox_cohort, hpo)

        phenotype_categories = (
            'HP:0012638',  # Abnormal nervous system physiology
            'HP:0001939',  # Abnormality of metabolism/homeostasis
        )
        variant_predicate = VariantPredicates.variant_effect(VariantEffect.MISSENSE_VARIANT, 'NM_001032386.2')

        counts_df = cohort_analysis.compare_symptom_count_vs_genotype(
            query=phenotype_categories, 
            variant_predicate=variant_predicate,
        )

        print(counts_df)
