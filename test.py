import unittest
from random import random
from typing import Callable

from math import factorial

from scipy.stats import fisher_exact

from fast_fisher import fast_fisher_python
from fast_fisher import fast_fisher_numba

try:
    from fast_fisher import fast_fisher_compiled
except ImportError:
    fast_fisher_numba.cc.compile()
    from fast_fisher import fast_fisher_compiled


def fisher_exact_factorial(a: int, b: int, c: int, d: int) -> float:
    """
    ((a+b)!*(c+d)!*(a+c)!*(b+d)!)
    -----------------------------
          (a!*b!*c!*d!*n!)
    """
    return (factorial(a + b) * factorial(c + d) * factorial(a + c) * factorial(b + d)) / \
           (factorial(a) * factorial(b) * factorial(c) * factorial(d) * factorial(a + b + c + d))


class TestFisher(unittest.TestCase):
    def test_with_scipy(self, samples=1000):
        def _test(function: Callable, name: str, table):
            l, r, t = function(*table)

            for pval in (l, r, t):
                self.assertLessEqual(pval, 1, f'Error in {name}: {pval=} must be smaller than 1! {table=}')

            self.assertAlmostEqual(l, scipy_l, msg=f'Error in {name}: {l=} != {scipy_l}; {table=}')
            self.assertAlmostEqual(r, scipy_r, msg=f'Error in {name}: {r=} != {scipy_r}; {table=}')
            self.assertAlmostEqual(t, scipy_t, msg=f'Error in {name}: {t=} != {scipy_t}; {table=}')

        for _ in range(samples):
            a, b, c, d = [int(10 ** (random() * 4) - 1) for _ in range(4)]

            scipy_l = fisher_exact([[a, b], [c, d]], 'less')[1]
            scipy_r = fisher_exact([[a, b], [c, d]], 'greater')[1]
            scipy_t = fisher_exact([[a, b], [c, d]], 'two-sided')[1]

            _test(fast_fisher_python.test1, 'fast_fisher_python', (a, b, c, d))
            _test(fast_fisher_numba.test1, 'fast_fisher_numba', (a, b, c, d))
            _test(fast_fisher_compiled.test1, 'fast_fisher_compiled', (a, b, c, d))

    def test_exception(self):
        a, b, c, d = (2, 3, 0, 2)
        scipy_t = fisher_exact([[a, b], [c, d]], 'two-sided')[1]
        fast_t = fast_fisher_python.test1t(a, b, c, d)

        # For this contingency table, scipy and this implementation do not agree!
        self.assertNotAlmostEqual(fast_t, scipy_t)

        factorial_t = fisher_exact_factorial(a, b, c, d)

        scipy_error = abs(scipy_t - factorial_t)
        fast_error = abs(fast_t - factorial_t)

        # This implementation gives a better result than scipy's implementation!
        self.assertTrue(fast_error < scipy_error)


if __name__ == '__main__':
    unittest.main()
