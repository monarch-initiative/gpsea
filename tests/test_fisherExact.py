import pytest
import numpy as np
from contextlib import nullcontext as does_not_raise
from genophenocorr.cohort._analyzers import PythonMultiFisherExact


@pytest.fixture
def MultiExact() -> PythonMultiFisherExact:
    return PythonMultiFisherExact()

@pytest.mark.parametrize('table, raise_error, pVal', 
                        ([[[0,0], [0,0], [0,0]], pytest.raises(ValueError), None],
                        [[[2,3], [1,0], [0,2]], does_not_raise(), 0.6429],
                        [[[100, 150], [500, 460], [420, 400]], pytest.raises(OverflowError), None],
                        [[[5,5],[5,5],[5,5]], does_not_raise(), 1],
                        [[[10,15], [5], [20,5]], pytest.raises(ValueError), None],
                        [[[10, 1], [2,3], [3,4]], does_not_raise(), 0.0395],
                        [[[6,6], [0,5], [0,8]], does_not_raise(), 0.0233],
                        [[[],[],[]], pytest.raises(ValueError), None]
))
def test_multiFisherExact(table, raise_error, pVal, MultiExact):
    with raise_error:
        np_table = np.array(table, dtype=np.int64)
        final_pval = MultiExact.calculate(np_table)
        assert round(final_pval, 4) == pVal
    