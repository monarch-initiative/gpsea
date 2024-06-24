import typing

import hpotk
import numpy as np
import pandas as pd
import pytest

from genophenocorr.analysis import HeuristicSamplerMtcFilter, apply_predicates_on_patients
from genophenocorr.analysis.predicate import PatientCategories, GenotypePolyPredicate
from genophenocorr.analysis.predicate.phenotype import PhenotypePolyPredicate
from genophenocorr.model import Cohort


class TestHeuristicSamplerMtcFilter:

    @pytest.fixture
    def mtc_filter(self, hpo: hpotk.MinimalOntology) -> HeuristicSamplerMtcFilter:
        return HeuristicSamplerMtcFilter(hpo=hpo)

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

        actual = HeuristicSamplerMtcFilter.one_genotype_has_zero_hpo_observations(counts=counts_df)

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

        actual = HeuristicSamplerMtcFilter.some_cell_has_greater_than_one_count(counts=counts_df)

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

        actual = HeuristicSamplerMtcFilter.genotypes_have_same_hpo_proportions(counts=counts_df)

        assert actual == expected

    def test_filter_terms_to_test(
            self,
            mtc_filter: HeuristicSamplerMtcFilter,
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
        assert reason_for_filtering_out['Skipping non-target term'] == 14
        assert reason_for_filtering_out['Skipping top level term'] == 5

        assert len(filtered_n_usable) == 4
        assert len(filtered_all_counts) == 4

    @staticmethod
    def prepare_counts_df(counts):
        index = pd.Index([PatientCategories.YES, PatientCategories.NO])
        columns = pd.Index([PatientCategories.YES, PatientCategories.NO])
        values = np.array(counts).reshape((2, 2))

        return pd.DataFrame(data=values, index=index, columns=columns)
