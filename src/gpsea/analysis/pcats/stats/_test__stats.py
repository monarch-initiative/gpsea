import numpy as np
import pandas as pd
import pytest

from ._stats import FisherExactTest


class TestPythonMultiFisherExact:

    @pytest.fixture
    def fisher_exact(self) -> FisherExactTest:
        return FisherExactTest()

    @pytest.mark.parametrize(
        "counts, expected",
        (
            [[[2, 1, 0], [3, 0, 2]], 0.6428571428571429],
            [[[5, 5, 5], [5, 5, 5]], 1],
            [[[10, 2, 3], [1, 3, 4]], 0.03952977071835599],
            [[[6, 0, 0], [6, 5, 8]], 0.02334274421230943],
        ),
    )
    def test_compute_pval(
        self,
        counts,
        expected,
        fisher_exact: FisherExactTest,
    ):
        contingency_matrix = pd.DataFrame(np.array(counts, dtype=np.int64))

        final_pval = fisher_exact.compute_pval(contingency_matrix)
        assert final_pval == pytest.approx(expected)
