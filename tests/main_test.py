import pytest
import sys
import pandas as pd
from genophenocorr import *


@pytest.fixture
def goodTestCohort():
    allPatients = AllPatients('testSamples/*.json', 'hg37')
    return allPatients

def test_all_properties(goodTestCohort):
    assert len(goodTestCohort.all_patients) == 10
    assert len(goodTestCohort.all_diseases) == 1
    assert len(goodTestCohort.all_phenotypes) == 62
    assert len(goodTestCohort.all_variants) == 10
    assert len(goodTestCohort.all_proteins) == 1
    assert len(goodTestCohort.all_var_types) == 5

def test_lists(goodTestCohort):
    assert goodTestCohort.list_all_diseases() == [['OMIM:616900', 'Hypotonia, infantile, with psychomotor retardation and characteristic facies 3']]
    assert goodTestCohort.list_all_proteins() == [['ENSP00000273980', 'TBC domain-containing protein kinase-like protein']]
    assert goodTestCohort.list_all_variants() == ['chr4 g.107115874C>T', 'chr4 g.107092429T>C', 'chr4 g.107156505_107156505delT', 'chr4 g.107152924A>G', 'chr4 g.107183260G>A', 'chr4 g.107183260G>A', 'chr4 g.107156505_107156505delT', 'chr4 g.107115874C>T', 'chr4 g.107183260G>A', 'chr4 g.107092429T>C']
    assert goodTestCohort.list_all_patients() == ['1-1', '4-2', '6-2', '3-1', '5-1', '2-1', '6-1', '1-2', '8-1', '4-1']

def test_count_per_feature(goodTestCohort):
    test_df = goodTestCohort.count_vars_per_feature()
    assert isinstance(test_df, pd.DataFrame)
    assert test_df.loc['ENSP00000273980'].at['Domain: Rab-GAP TBC', 'variants'] == 1
