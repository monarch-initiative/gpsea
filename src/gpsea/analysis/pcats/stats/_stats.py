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

    
    Supports shape
    ^^^^^^^^^^^^^^

    `CountStatistic` takes the counts in form of a data frame,
    and some statistics impose additional requirements on the frame shape.
    For instance, GPSEA's implementation of the Fisher exact test
    can compare counts in a ``(2, 2)`` or ``(2, 3)`` arrays
    but χ2 test can test an ``(m, n)`` array.

    It is important to check that a genotype/phenotype predicate produces
    the number of groups which the statistic can test.

    The :attr:`supports_shape` returns a sequence with requirements
    on the shape of the data array/frame. The sequence includes
    the number of

    Examples
    ********

    +------------------------+-------------------------+------------------+
    | Test                   | Array shape             | `supports_shape` |
    +========================+=========================+==================+
    | Fisher Exact Test      | ``(2, [2, 3])``         | ``(2, [2,3])``   |
    +------------------------+-------------------------+------------------+
    | χ2                     | ``(*, *)``              | ``(None, None)`` |
    +------------------------+-------------------------+------------------+
    """

    @property
    @abc.abstractmethod
    def supports_shape(
        self,
    ) -> typing.Sequence[typing.Union[int, typing.Sequence[int], None]]:
        """
        Get a sequence of the supported shapes.
        """
        pass

    @abc.abstractmethod
    def compute_pval(
        self,
        counts: pd.DataFrame,
    ) -> float:
        pass


class FisherExactTest(CountStatistic):
    """
    `FisherExactTest` performs Fisher Exact Test on a `2x2` or `2x3` contingency table.

    The `2x2` version is a thin wrapper around Scipy :func:`~scipy.stats.fisher_exact` function,
    while the `2x3` variant is implemented in Python.
    In both variants, the two-sided :math:`H_1` is considered.
    """
    
    def __init__(self):
        self._shape = (2, (2, 3))

    @property
    def supports_shape(
        self,
    ) -> typing.Sequence[typing.Union[int, typing.Sequence[int], None]]:
        return self._shape

    def compute_pval(
        self,
        counts: pd.DataFrame,
    ) -> float:
        if counts.shape == (2, 2):
            _, pval = scipy.stats.fisher_exact(counts.values, alternative="two-sided")
            return pval
        elif counts.shape == (2, 3):
            return self._fisher_exact(counts.values)
        else:
            raise ValueError(f'Unsupported counts shape {counts.shape}')

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
