import abc
import math
import typing
from decimal import Decimal

import numpy as np
import scipy


class MultiFisherExact(metaclass=abc.ABCMeta):

    @abc.abstractmethod
    def calculate(self, a: np.ndarray) -> float:
        """
        :param a: a 2x3 int array with counts
        :returns: a p value calculated with Fisher's exact test
        """
        pass

    @staticmethod
    def _check_input(a: np.ndarray):
        if not isinstance(a, np.ndarray):
            raise ValueError(f'Expected a numpy array but got {type(a)}')
        if not a.shape == (2, 3):
            raise ValueError(f'Shape of the array must be (2, 3) but got {a.shape}')
        if np.array_equal(a, np.zeros_like(a)):
            raise ValueError(f'Array is all zeros, cannot run analysis')


class PythonMultiFisherExact(MultiFisherExact):

    def calculate(self, a: np.ndarray) -> float:
        MultiFisherExact._check_input(a)
        return self._fisher_exact(a)

    def _fisher_exact(self, table):
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


def run_recessive_fisher_exact(two_by_three_table: typing.Sequence[typing.Sequence[int]]):
    a = np.array(two_by_three_table, dtype=np.int64)
    test_class = PythonMultiFisherExact()
    val = test_class.calculate(a)
    return val


def run_fisher_exact(two_by_two_table: typing.Sequence[typing.Sequence[int]]):
    oddsr, p = scipy.stats.fisher_exact(two_by_two_table, alternative='two-sided')
    return p
