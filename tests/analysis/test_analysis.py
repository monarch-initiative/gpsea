import typing

import hpotk
import numpy as np
import pandas as pd

from genophenocorr.analysis import apply_predicates_on_patients

from genophenocorr.model import Cohort
from genophenocorr.analysis.predicate import GenotypePolyPredicate, PatientCategories
from genophenocorr.analysis.predicate.phenotype import PhenotypePolyPredicate


def test_apply_predicates_on_patients(
        suox_cohort: Cohort,
        suox_pheno_predicates: typing.Sequence[PhenotypePolyPredicate[hpotk.TermId]],
        suox_gt_predicate: GenotypePolyPredicate,
):
    categories, n_usable, counts = apply_predicates_on_patients(
        patients=suox_cohort.all_patients,
        pheno_predicates=suox_pheno_predicates,
        gt_predicate=suox_gt_predicate,
    )

    assert isinstance(categories, typing.Collection)
    assert PatientCategories.YES in categories
    assert PatientCategories.NO in categories
    assert len(categories) == 2

    assert isinstance(n_usable, pd.Series)
    assert len(n_usable) == 5

    assert isinstance(counts, typing.Mapping)
    assert len(counts) == 5

    seizure_counts = counts[hpotk.TermId.from_curie('HP:0001250')]
    assert np.array_equal(seizure_counts.to_numpy(), np.array([[17, 11], [7, 0]]))
