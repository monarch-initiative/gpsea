import pytest
import pandas as pd

@pytest.fixture
def outputVars():
    out = pd.read_csv('testSamples/compare_outputs.tsv', sep='\t')

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
def test_compare(patient, num):
    pass
