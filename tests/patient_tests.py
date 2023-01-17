import pytest
import sys
import pandas as pd
from genophenocorr import *

@pytest.fixture
def outputVars():
    out = pd.read_csv('testSamples/testing_outputs.tsv', sep='\t')
    out['disease'] = out['disease'].apply(eval)
    out['protein'] = out['protein'].apply(eval)
    out['phenotypes'] = out['phenotypes'].apply(eval)
    return out

@pytest.mark.parametrize('patient, num',
    [('testSamples/Bhoj-2016-TBCK-1-1.json', 0),
    ('testSamples/Bhoj-2016-TBCK-1-2.json', 1),
    ('testSamples/Bhoj-2016-TBCK-2-1.json', 2),
    ('testSamples/Bhoj-2016-TBCK-3-1.json', 3),
    ('testSamples/Bhoj-2016-TBCK-4-1.json', 4),
    ('testSamples/Bhoj-2016-TBCK-4-2.json', 5),
    ('testSamples/Bhoj-2016-TBCK-5-1.json', 6),
    ('testSamples/Bhoj-2016-TBCK-6-1.json', 7),
    ('testSamples/Bhoj-2016-TBCK-6-2.json', 8),
    ('testSamples/Bhoj-2016-TBCK-8-1.json', 9)
    ])

def test_patient(patient, num, outputVars):
    patTest = Patient(patient, 'hg37')
    assert isinstance(patTest, Patient)
    assert patTest.id == outputVars.at[num, 'id']
    assert patTest.variant.variant_string == outputVars.at[num, 'variant']
    assert patTest.disease_id == outputVars.at[num, 'disease'][0]
    assert patTest.disease_label == outputVars.at[num, 'disease'][1]
    assert patTest.protein.id == outputVars.at[num, 'protein'][0]
    assert patTest.protein.label == outputVars.at[num, 'protein'][1]
    assert patTest.phenotype_ids == outputVars.at[num, 'phenotypes']