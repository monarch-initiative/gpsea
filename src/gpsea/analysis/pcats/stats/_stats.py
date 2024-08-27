import abc
import math
import typing
from decimal import Decimal


import numpy as np
import pandas as pd

import scipy


class CountStatistic(metaclass=abc.ABCMeta):
    """
    `CountStatistic` calculates a p value for a contingency table
    produced by a pair of discrete random variables.

    The `counts` table is usually `2x2` or `2x3`.
    """

    @abc.abstractmethod
    def compute_pval(
        self,
        counts: pd.DataFrame,
    ) -> float:
        pass


class ScipyFisherExact(CountStatistic):
    """
    `ScipyFisherExact` performs Fisher Exact Test on a `2x2` contingency table.

    The class is a thin wrapper around Scipy :func:`~scipy.stats.fisher_exact` function.
    The two-sided :math:`H_1` is considered.
    """

    def compute_pval(
        self,
        counts: pd.DataFrame,
    ) -> float:
        assert counts.shape == (
            2,
            2,
        ), f"Cannot run Fisher exact test on an array with {counts.shape} shape"
        _, pval = scipy.stats.fisher_exact(counts.values, alternative="two-sided")
        return pval


class PythonMultiFisherExact(CountStatistic):
    """
    `PythonMultiFisherExact` is a Python implementation of Fisher Exact Test to compute
    p value for a `2x3` contingency table.
    """

    def compute_pval(
        self,
        counts: pd.DataFrame,
    ) -> float:
        PythonMultiFisherExact._check_input(counts)
        return self._fisher_exact(counts.values)

    @staticmethod
    def _check_input(a: pd.DataFrame):
        if not isinstance(a, pd.DataFrame):
            raise ValueError(f"Expected a pandas DataFrame but got {type(a)}")
        if not a.shape == (2, 3):
            raise ValueError(f"Shape of the array must be (2, 3) but got {a.shape}")
        if np.array_equal(a.values, np.zeros_like(a)):
            raise ValueError("Data frame is all zeros, cannot run analysis")

    def _fisher_exact(
        self,
        table: np.ndarray,
    ):
        row_sum = []
        col_sum = []

        for i in range(len(table)):
            temp = 0
            for j in range(len(table[0])):
                temp += table[i][j]
            row_sum.append(temp)

        for j in range(len(table[0])):
            temp = 0
            for i in range(len(table)):
                temp += table[i][j]
            col_sum.append(temp)

        mat = [[0] * len(col_sum)] * len(row_sum)
        pos = (0, 0)

        p_0 = 1

        for x in row_sum:
            p_0 *= math.factorial(x)
        for y in col_sum:
            p_0 *= math.factorial(y)

        n = 0
        for x in row_sum:
            n += x
        p_0 /= Decimal(math.factorial(n))

        for i in range(len(table)):
            for j in range(len(table[0])):
                p_0 /= Decimal(math.factorial(table[i][j]))

        p = [0]
        self._dfs(mat, pos, row_sum, col_sum, p_0, p)

        return float(p[0])

    def _dfs(self, mat, pos, r_sum, c_sum, p_0, p):

        (xx, yy) = pos
        (r, c) = (len(r_sum), len(c_sum))

        mat_new = []

        for i in range(len(mat)):
            temp = []
            for j in range(len(mat[0])):
                temp.append(mat[i][j])
            mat_new.append(temp)

        if xx == -1 and yy == -1:
            for i in range(r - 1):
                temp = r_sum[i]
                for j in range(c - 1):
                    temp -= mat_new[i][j]
                mat_new[i][c - 1] = temp
            for j in range(c - 1):
                temp = c_sum[j]
                for i in range(r - 1):
                    temp -= mat_new[i][j]
                mat_new[r - 1][j] = temp
            temp = r_sum[r - 1]
            for j in range(c - 1):
                temp -= mat_new[r - 1][j]
            if temp < 0:
                return
            mat_new[r - 1][c - 1] = temp

            p_1 = 1
            for x in r_sum:
                p_1 *= math.factorial(x)
            for y in c_sum:
                p_1 *= math.factorial(y)

            n = 0
            for x in r_sum:
                n += x
            p_1 /= Decimal(math.factorial(n))

            for i in range(len(mat_new)):
                for j in range(len(mat_new[0])):
                    p_1 /= Decimal(math.factorial(mat_new[i][j]))
            if p_1 <= p_0 + Decimal(0.00000001):
                # print(mat_new)
                # print(p_1)
                p[0] += p_1
        else:
            max_1 = r_sum[xx]
            max_2 = c_sum[yy]
            for j in range(c):
                max_1 -= mat_new[xx][j]
            for i in range(r):
                max_2 -= mat_new[i][yy]
            for k in range(min(max_1, max_2) + 1):
                mat_new[xx][yy] = k
                if xx == r - 2 and yy == c - 2:
                    pos_new = (-1, -1)
                elif xx == r - 2:
                    pos_new = (0, yy + 1)
                else:
                    pos_new = (xx + 1, yy)
                self._dfs(mat_new, pos_new, r_sum, c_sum, p_0, p)
