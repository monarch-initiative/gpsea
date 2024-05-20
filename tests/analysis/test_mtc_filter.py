import typing

import hpotk
import pytest

import pandas as pd
import numpy as np

from genophenocorr.analysis import HeuristicSamplerMtcFilter
from genophenocorr.analysis.predicate import PatientCategories

class TestHeuristicSamplerMtcFilter():
    
    @pytest.fixture
    def mtc_filter(self, hpo: hpotk.MinimalOntology) -> HeuristicSamplerMtcFilter:
        return HeuristicSamplerMtcFilter(hpo_ontology=hpo)
    
    
    @pytest.mark.parametrize(
        "counts, expected",
        [
            ((95,60, 144-95,71-60), False),
            ((40,15,18,15), False)
        ]
        )
    def test_genotypes_have_same_HPO_proportions(
        self, 
        counts: typing.Tuple[int], 
        expected: bool,
    ):
        index = pd.Index([PatientCategories.YES, PatientCategories.NO])
        columns = pd.Index([PatientCategories.YES, PatientCategories.NO])
        values = np.array(counts).reshape((2,2))
        counts_df = pd.DataFrame(data=values, index=index, columns=columns)
        actual = HeuristicSamplerMtcFilter.genotypes_have_same_HPO_proportions(counts=counts_df)
        assert actual == expected