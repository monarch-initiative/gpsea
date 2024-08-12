import typing

import hpotk
import numpy as np
import pandas as pd
import pytest

from genophenocorr.analysis import HeuristicMtcFilter, SpecifiedTermsMtcFilter, apply_predicates_on_patients
from genophenocorr.analysis.predicate import PatientCategories, GenotypePolyPredicate
from genophenocorr.analysis.predicate.phenotype import PhenotypePolyPredicate
from genophenocorr.model import Cohort


class TestHeuristicSamplerMtcFilter:

    @pytest.fixture
    def mtc_filter(self, hpo: hpotk.MinimalOntology) -> HeuristicMtcFilter:
        default_freq_threshold=0.2
        return HeuristicMtcFilter.default_filter(
            hpo=hpo, 
            term_frequency_threshold=default_freq_threshold,
        )

    @pytest.fixture(scope='class')
    def patient_counts(
            self,
            suox_cohort: Cohort,
            suox_gt_predicate: GenotypePolyPredicate,
            suox_pheno_predicates: typing.Sequence[PhenotypePolyPredicate[hpotk.TermId]],
    ) -> typing.Tuple[
        typing.Mapping[hpotk.TermId, int],
        typing.Mapping[hpotk.TermId, pd.DataFrame],
    ]:
        categories, n_usable, all_counts = apply_predicates_on_patients(
            patients=suox_cohort.all_patients,
            gt_predicate=suox_gt_predicate,
            pheno_predicates=suox_pheno_predicates,
        )
        return n_usable, all_counts

    @pytest.mark.parametrize(
        "counts, expected",
        [
            ((10, 1, 15, 0), False),
            ((10, 0, 15, 1), False),
            ((10, 0, 15, 0), True),
            ((0, 15, 0, 19), True),
        ]
    )
    def test_one_genotype_has_zero_hpo_observations(
            self,
            counts: typing.Tuple[int],
            expected: bool,
    ):
        counts_df = self.prepare_counts_df(counts)

        actual = HeuristicMtcFilter.one_genotype_has_zero_hpo_observations(counts=counts_df)

        assert actual == expected

    @pytest.mark.parametrize(
        "counts, expected",
        [
            ((1, 0, 20, 30), False),
            ((2, 0, 20, 30), True),

            ((0, 1, 20, 30), False),
            ((0, 2, 20, 30), True),

            ((0, 0, 20, 30), False),
            ((2, 2, 20, 30), True),
        ]
    )
    def test_some_cell_has_greater_than_one_count(
            self,
            counts: typing.Tuple[int],
            expected: bool,
    ):
        counts_df = self.prepare_counts_df(counts)

        actual = HeuristicMtcFilter.some_cell_has_greater_than_one_count(counts=counts_df)

        assert actual == expected

    @pytest.mark.parametrize(
        "counts, expected",
        [
            ((10, 100, 50, 500), True),
            ((95, 60, 144 - 95, 71 - 60), False),
            ((40, 15, 18, 15), False),
        ]
    )
    def test_genotypes_have_same_hpo_proportions(
            self,
            counts: typing.Tuple[int],
            expected: bool,
    ):
        counts_df = self.prepare_counts_df(counts)

        actual = HeuristicMtcFilter.genotypes_have_same_hpo_proportions(counts=counts_df)

        assert actual == expected

    def test_filter_terms_to_test(
            self,
            mtc_filter: HeuristicMtcFilter,
            patient_counts: typing.Tuple[
                typing.Mapping[hpotk.TermId, int],
                typing.Mapping[hpotk.TermId, pd.DataFrame],
            ],
    ):
        n_usable, all_counts = patient_counts
        mtc_report = mtc_filter.filter_terms_to_test(n_usable, all_counts)

        assert isinstance(mtc_report, tuple)
        assert len(mtc_report) == 3

        filtered_n_usable, filtered_all_counts, reason_for_filtering_out = mtc_report

        assert reason_for_filtering_out['Skipping term because all genotypes have same HPO observed proportions'] == 1
        assert reason_for_filtering_out['Skipping general term'] == 14
        assert reason_for_filtering_out['Skipping non-target term'] == 5
        assert reason_for_filtering_out['Skipping top level term'] == 0

        assert len(filtered_n_usable) == 4
        assert len(filtered_all_counts) == 4

    def test_specified_term_mtc_filter(
            self,
            hpo: hpotk.MinimalOntology,
            patient_counts: typing.Tuple[
                typing.Mapping[hpotk.TermId, int],
                typing.Mapping[hpotk.TermId, pd.DataFrame],
            ],
    ):
        """
        The point of this test is to check that if we filter to test only one term ("HP:0032350"), then this
        is the only term that should survive the filter. We start with a total of five terms (n_usable==5),
        but after our filter, only one survives (filtered_n_usable == 1), and we have four cases in which the
        reason for filtering out is 'Skipping non-specified term'
        """
        specified_filter = SpecifiedTermsMtcFilter(hpo=hpo, terms_to_test={hpotk.TermId.from_curie("HP:0032350")})
        n_usable, all_counts = patient_counts
        mtc_report = specified_filter.filter_terms_to_test(n_usable, all_counts)
        assert isinstance(mtc_report, tuple)
        assert len(mtc_report) == 3 ## Skipping non-specified term (n=5)

        filtered_n_usable, filtered_all_counts, reason_for_filtering_out = mtc_report
        assert len(n_usable) == 5
        assert len(filtered_n_usable) == 1
        assert reason_for_filtering_out['Skipping non-specified term'] == 4


    @staticmethod
    def prepare_counts_df(counts):
        index = pd.Index([PatientCategories.YES, PatientCategories.NO])
        columns = pd.Index([PatientCategories.YES, PatientCategories.NO])
        values = np.array(counts).reshape((2, 2))

        return pd.DataFrame(data=values, index=index, columns=columns)


    def test_min_observed_HPO_threshold(
        self,
        patient_counts: typing.Tuple[
            typing.Mapping[hpotk.TermId, int],
            typing.Mapping[hpotk.TermId, pd.DataFrame],
            ],
    ):
        """
        In our heuristic filter, we only test terms that have at least a threshold
        frequency in at least one of the groups. We use the "all counts" datastructure, that
        is a dictionary whose keys are hpotoolkit TermIds and whose values are pandas DataFrames
        with 2x2 contingenicy tables of counts. For instance, each column will have one row for
        PatientCategories.YES and one for PatientCategories.NO, indicating counts of measured observed/excluded
        HPO phenotypes. Each column is a certain genotype, e.g., MISSENSE or NON-MISSENSE. We want the
        function to return the maximum frequency. In each column, the frequency is calculate by
        PatientCategories.YES / (PatientCategories.YES+PatientCategories.NO). This function tests that this works
        for all of the HPO terms in the dictionary.
        """
        EPSILON = 0.001
        _, all_counts = patient_counts
        # Ectopia lentis HP:0001083  (6 9  1 2), freqs are 6/15=0.4 and 1/3=0.33
        ectopia = all_counts[hpotk.TermId.from_curie("HP:0001083")]
        max_f = HeuristicMtcFilter.get_maximum_group_observed_HPO_frequency(ectopia)
        assert max_f == pytest.approx(0.4, abs=EPSILON)
        
        # Seizure HP:0001250 (17 7 11 0), freqs are 17/24=0.7083 and 11/11=1
        seizure = all_counts[hpotk.TermId.from_curie("HP:0001250")]
        max_f = HeuristicMtcFilter.get_maximum_group_observed_HPO_frequency(seizure)
        assert max_f == pytest.approx(1.0, abs=EPSILON)
        
        # Sulfocysteinuria HP:0032350 (11 0 2 0), freqs are both 1
        sulfocysteinuria = all_counts[hpotk.TermId.from_curie("HP:0032350")]
        max_f = HeuristicMtcFilter.get_maximum_group_observed_HPO_frequency(sulfocysteinuria)
        assert max_f == pytest.approx(1.0, abs=EPSILON)
        
        # Neurodevelopmental delay HP:0012758 (4 13 4 4), freqs are 4/17 = 0.235 and 4/8=0.5
        ndelay = all_counts[hpotk.TermId.from_curie("HP:0012758")]
        max_f = HeuristicMtcFilter.get_maximum_group_observed_HPO_frequency(ndelay)
        assert max_f == pytest.approx(0.5, abs=EPSILON)
        
        # Hypertonia HP:0001276 (7 9 4 3) fresa are 7/16=0.4375 and 4/7=0.5714
        hypertonia  = all_counts[hpotk.TermId.from_curie("HP:0001276")]
        max_f = HeuristicMtcFilter.get_maximum_group_observed_HPO_frequency(hypertonia)
        assert max_f == pytest.approx(0.5714, abs=EPSILON)



