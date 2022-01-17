#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: eph

from timeit import timeit
from scipy import stats
from fast_fisher import fast_fisher, fast_fisher_numba

if __name__ == '__main__':

    test_mapping = {'left-tailed': 'less',
                    'right-tailed': 'greater',
                    'two-tailed': 'two-sided'}

    setup = 'from __main__ import stats, fast_fisher, fast_fisher_numba'

    print('| {:>6s} | {:>6s} | {:>6s} | {:>6s} | {:>12s} | {:>9s} | {:>9s} | {:>9s} |'
          .format('a', 'b', 'c', 'd', 'test type', 'scipy', 'fisher', 'f_numba'))
    print('|-------:|-------:|-------:|-------:|-------------:|----------:|----------:|----------:|')
    for a, b, c, d, test_type in [
        (8, 2, 1, 5, 'left-tailed'),
        (8, 2, 1, 5, 'right-tailed'),
        (8, 2, 1, 5, 'two-tailed'),
        (8, 2, 1, 5, 'left-tailed'),
        (8, 2, 1, 5, 'right-tailed'),
        (8, 2, 1, 5, 'two-tailed'),
        (100, 1000, 10000, 100000, 'left-tailed'),
        (100, 1000, 10000, 100000, 'right-tailed'),
        (100, 1000, 10000, 100000, 'two-tailed'),
        (10000, 100, 1000, 100000, 'left-tailed'),
        (10000, 100, 1000, 100000, 'right-tailed'),
        (10000, 100, 1000, 100000, 'two-tailed'),
        (10000, 10000, 10000, 10000, 'left-tailed'),
        (10000, 10000, 10000, 10000, 'right-tailed'),
        (10000, 10000, 10000, 10000, 'two-tailed'),
    ]:
        print('| {:>6d} | {:>6d} | {:>6d} | {:>6d} | {:>12s} | {:>6.0f} us | {:>6.0f} us | {:>6.0f} us |'.format(
            a, b, c, d, test_type,
            timeit('stats.fisher_exact([[{}, {}], [{}, {}]], "{}")'
                   .format(a, b, c, d, test_mapping[test_type]),
                   setup=setup, number=100) * 1e4,
            timeit('fast_fisher.test1{}({}, {}, {}, {})'
                   .format(test_type[0], a, b, c, d),
                   setup=setup, number=100) * 1e4,
            timeit('fast_fisher_numba.test1{}({}, {}, {}, {})'
                   .format(test_type[0], a, b, c, d),
                   setup=setup, number=100) * 1e4,
        ))
