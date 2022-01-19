__version__ = '0.0.2'

from . import fast_fisher_python

try:
    # try to import compiled version
    from . import fast_fisher_compiled

    fast_fisher = fast_fisher_compiled
except ImportError as e:
    # import python version
    from logging import warning

    warning(f'Failed to import fast_fisher_compiled: {str(e)}')
    fast_fisher = fast_fisher_python

odds_ratio = fast_fisher.odds_ratio


def fast_fisher_exact(a: int, b: int, c: int, d: int, alternative: str = 'two-sided'):
    """
    Perform a Fisher exact test on a 2x2 contingency table.

    :param a: row 1 col 1
    :param b: row 1 col 2
    :param c: row 2 col 1
    :param d: row 2 col 2
    :param alternative: {‘two-sided’, ‘less’, ‘greater’} (default: 'two-sided')
    :return: pvalue
    """
    if alternative is None:
        alternative = 'two-sided'

    return fast_fisher.fisher_exact(a, b, c, d, alternative)


def fast_fisher_exact_compatibility(table: [[int, int], [int, int]], alternative: str = 'two-sided'):
    """
    Perform a Fisher exact test on a 2x2 contingency table.

    :param table: A 2x2 contingency table. Elements must be non-negative integers.
    :param alternative: {‘two-sided’, ‘less’, ‘greater’} (default: 'two-sided')
    :return: pvalue
    """
    (a, b), (c, d) = table
    return odds_ratio(a, b, c, d), fast_fisher.fisher_exact(a, b, c, d, alternative)


try:
    from . import fast_fisher_numba
except ImportError as e:
    from logging import warning

    warning(f'Failed to import fast_fisher_numba: {str(e)}')
    pass
