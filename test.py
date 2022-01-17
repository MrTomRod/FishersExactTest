#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: eph

import unittest
from random import random

from fast_fisher import fast_fisher
from fast_fisher import fast_fisher_numba


class TestFisher(unittest.TestCase):
    def test_with_scipy(self, samples=100):
        from scipy.stats import fisher_exact
        for _ in range(samples):
            a, b, c, d = [int(10 ** (random() * 4) - 1) for _ in range(4)]

            scipy_l = fisher_exact([[a, b], [c, d]], 'less')[1]
            scipy_r = fisher_exact([[a, b], [c, d]], 'greater')[1]
            scipy_t = fisher_exact([[a, b], [c, d]], 'two-sided')[1]

            fast_l, fast_r, fast_t = fast_fisher.test1(a, b, c, d)
            numba_l, numba_r, numba_t = fast_fisher_numba.test1(a, b, c, d)

            for pval in (fast_l, fast_r, fast_t, numba_l, numba_r, numba_t):
                self.assertLessEqual(pval, 1)

            for test1 in (fast_fisher.test1, fast_fisher_numba.test1):
                l, r, t = test1(a, b, c, d)
                self.assertAlmostEqual(l, scipy_l)
                self.assertAlmostEqual(r, scipy_r)
                self.assertAlmostEqual(t, scipy_t)


if __name__ == '__main__':
    unittest.main()
