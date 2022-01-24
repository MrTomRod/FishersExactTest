from unittest import TestCase

from random import random
from typing import Callable
from itertools import combinations, product

from math import isclose, isnan, isinf
from random import randint

import pandas as pd

from scipy.stats import fisher_exact

from fast_fisher import fast_fisher_exact, fast_fisher_exact_compatibility, fast_fisher_python, fast_fisher_numba, fast_fisher_cython

try:
    from fast_fisher import fast_fisher_compiled
except ImportError:
    fast_fisher_numba.cc.compile()
    from fast_fisher import fast_fisher_compiled

FISHER_COMBINATIONS = ['abcd', 'acbd', 'badc', 'bdac', 'cadb', 'cdab', 'dbca', 'dcba']


def scipy_fisher_exact(a, b, c, d, alternative: str = 'two-sided'):
    return fisher_exact([[a, b], [c, d]], alternative=alternative)[1]


def fisher_swap(a: int, b: int, c: int, d: int) -> (int, int, int, int):
    """
    Eight contingency tables always give the same pvalue: ['abcd', 'acbd', 'badc', 'bdac', 'cadb', 'cdab', 'dbca', 'dcba']

    Compute and save only one version.
    """
    vals = {'a': a, 'b': b, 'c': c, 'd': d}
    equivalent_combinations = [tuple(vals[letter] for letter in combination) for combination in FISHER_COMBINATIONS]
    return sorted(equivalent_combinations, key=lambda comb: (comb[0], comb[1], comb[2], comb[3]), reverse=True)[0]


def get_disagreement_df(res):
    df = pd.DataFrame(res).T
    df.columns = ['scipy', 'fast']
    df = df[df.scipy + df.fast < 1.999999]
    df['rank_scipy'] = df.scipy.rank(method='min')
    df['rank_fast'] = df.fast.rank(method='min')
    df['agreement'] = df.rank_scipy == df.rank_fast
    df.sort_values(by='rank_scipy', inplace=True)

    disagreement_df = df[df.rank_scipy != df.rank_fast]
    first_disagreement_rank = disagreement_df.iloc[0]
    return disagreement_df, first_disagreement_rank


def is_equivalent(a: float, b: float) -> bool:
    if isinf(a) and isinf(b):
        return True
    if isnan(a) and isnan(b):
        return True
    return isclose(a, b)
