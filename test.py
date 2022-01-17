import unittest
from random import random
from typing import Callable

from scipy.stats import fisher_exact

from fast_fisher import fast_fisher_python
from fast_fisher import fast_fisher_numba

try:
    from fast_fisher import fast_fisher_compiled
except ImportError:
    fast_fisher_numba.cc.compile()
    from fast_fisher import fast_fisher_compiled


class TestFisher(unittest.TestCase):
    def _test(self, function: Callable, name: str, args, scipy_l, scipy_r, scipy_t):
        l, r, t = function(*args)

        for pval in (l, r, t):
            self.assertLessEqual(pval, 1, f'Error in {name}: {pval=} must be smaller than 1! {args=}')

        self.assertAlmostEqual(l, scipy_l, msg=f'Error in {name}: {l} != {scipy_l}; {args=}')
        self.assertAlmostEqual(r, scipy_r, msg=f'Error in {name}: {r} != {scipy_r}; {args=}')
        self.assertAlmostEqual(t, scipy_t, msg=f'Error in {name}: {l} != {scipy_t}; {args=}')

    def test_with_scipy(self, samples=100):
        for _ in range(samples):
            a, b, c, d = [int(10 ** (random() * 4) - 1) for _ in range(4)]

            scipy_l = fisher_exact([[a, b], [c, d]], 'less')[1]
            scipy_r = fisher_exact([[a, b], [c, d]], 'greater')[1]
            scipy_t = fisher_exact([[a, b], [c, d]], 'two-sided')[1]

            self._test(fast_fisher_python.test1, 'fast_fisher_python', (a, b, c, d), scipy_l, scipy_r, scipy_t)
            self._test(fast_fisher_numba.test1, 'fast_fisher_numba', (a, b, c, d), scipy_l, scipy_r, scipy_t)
            self._test(fast_fisher_compiled.test1, 'fast_fisher_compiled', (a, b, c, d), scipy_l, scipy_r, scipy_t)


if __name__ == '__main__':
    unittest.main()
