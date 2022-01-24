# cython: boundscheck=False
# cython: cdivision=True
# cython: wraparound=False
# cython: infer_types=True

from math import nan as pymath_nan, inf as pymath_inf
from libc.math cimport log, exp, lgamma, INFINITY, llround

cdef inline _maxn():
    l, n, h = 1, 2, INFINITY
    while l < n:
        if abs(lgamma(n + 1) - lgamma(n) - log(n)) >= 1:
            h = n
        else:
            l = n
        n = (l + min(h, l * 3)) // 2
    return n

cdef double LN10 = log(10)
cdef double NINF = INFINITY
cdef  long long MAXN = _maxn()

# ======================== Full Test ========================
cpdef (double, double, double) test1(long long a, long long b, long long c, long long d) except *:
    result = mlnTest2(a, a + b, a + c, a + b + c + d)
    return exp(-result[0]), exp(-result[1]), exp(-result[2])

cpdef (double, double, double) test2(long long a, long long ab, long long ac, long long abcd) except *:
    result = mlnTest2(a, ab, ac, abcd)
    return exp(-result[0]), exp(-result[1]), exp(-result[2])

cdef inline (double, double, double) mlnTest1(long long a, long long b, long long c, long long d) except *:
    return mlnTest2(a, a + b, a + c, a + b + c + d)

cdef inline (double, double, double) mlnTest2(long long a, long long ab, long long ac, long long abcd) except *:
    if 0 > a or a > ab or a > ac or ab + ac > abcd + a:
        raise ValueError('invalid contingency table')
    if abcd > MAXN:
        raise OverflowError('the grand total of contingency table is too large')
    cdef long long a_min = max(0, ab + ac - abcd)
    cdef long long a_max = min(ab, ac)
    if a_min == a_max:
        return 0., 0., 0.
    cdef double p0 = lgamma(ab + 1) + lgamma(ac + 1) + lgamma(abcd - ac + 1) + lgamma(abcd - ab + 1) - lgamma(abcd + 1)
    cdef double pa = lgamma(a + 1) + lgamma(ab - a + 1) + lgamma(ac - a + 1) + lgamma(abcd - ab - ac + a + 1)
    cdef double sl = 0.
    cdef double sr = 0.
    cdef double pi, sl_new, sr_new
    if ab * ac < a * abcd:
        for i in range(min(a - 1, llround(ab * ac / abcd)), a_min - 1, -1):
            pi = lgamma(i + 1) + lgamma(ab - i + 1) + lgamma(ac - i + 1) + lgamma(abcd - ab - ac + i + 1)
            if pi < pa:
                continue
            sl_new = sl + exp(pa - pi)
            if sl_new == sl:
                break
            sl = sl_new
        for i in range(a + 1, a_max + 1):
            pi = lgamma(i + 1) + lgamma(ab - i + 1) + lgamma(ac - i + 1) + lgamma(abcd - ab - ac + i + 1)
            sr_new = sr + exp(pa - pi)
            if sr_new == sr:
                break
            sr = sr_new
        return -log(1. - max(0, exp(p0 - pa) * sr)), max(0, pa - p0 - log(1. + sr)), max(0, pa - p0 - log(sl + 1. + sr))
    else:
        for i in range(a - 1, a_min - 1, -1):
            pi = lgamma(i + 1) + lgamma(ab - i + 1) + lgamma(ac - i + 1) + lgamma(abcd - ab - ac + i + 1)
            sl_new = sl + exp(pa - pi)
            if sl_new == sl:
                break
            sl = sl_new
        for i in range(max(a + 1, llround(ab * ac / abcd)), a_max + 1):
            pi = lgamma(i + 1) + lgamma(ab - i + 1) + lgamma(ac - i + 1) + lgamma(abcd - ab - ac + i + 1)
            if pi < pa:
                continue
            sr_new = sr + exp(pa - pi)
            if sr_new == sr:
                break
            sr = sr_new
        return max(0, pa - p0 - log(sl + 1.)), -log(1. - max(0, exp(p0 - pa) * sl)), max(0, pa - p0 - log(sl + 1. + sr))

cdef inline (double, double, double) mlog10Test1(long long a, long long b, long long c, long long d) except *:
    cdef double r1, r2, r3
    r1, r2, r3 = mlnTest2(a, a + b, a + c, a + b + c + d)
    return r1 / LN10, r2 / LN10, r3 / LN10

cdef inline (double, double, double) mlog10Test2(long long a, long long ab, long long ac, long long abcd) except *:
    result = mlnTest2(a, ab, ac, abcd)
    return result[0] / LN10, result[1] / LN10, result[2] / LN10

# ======================== Left Tail Only ========================
cpdef double test1l(long long a, long long b, long long c, long long d) except *:
    return exp(-mlnTest2l(a, a + b, a + c, a + b + c + d))

cpdef double test2l(long long a, long long ab, long long ac, long long abcd) except *:
    return exp(-mlnTest2l(a, ab, ac, abcd))

cdef inline double mlnTest1l(long long a, long long b, long long c, long long d) except *:
    return mlnTest2l(a, a + b, a + c, a + b + c + d)

cdef inline double mlnTest2l(long long a, long long ab, long long ac, long long abcd) except *:
    if 0 > a or a > ab or a > ac or ab + ac > abcd + a:
        raise ValueError('invalid contingency table')
    if abcd > MAXN:
        raise OverflowError('the grand total of contingency table is too large')
    cdef long long a_min = max(0, ab + ac - abcd)
    cdef long long a_max = min(ab, ac)
    if a_min == a_max:
        return 0.
    cdef double p0 = lgamma(ab + 1) + lgamma(ac + 1) + lgamma(abcd - ac + 1) + lgamma(abcd - ab + 1) - lgamma(abcd + 1)
    cdef double pa = lgamma(a + 1) + lgamma(ab - a + 1) + lgamma(ac - a + 1) + lgamma(abcd - ab - ac + a + 1)
    cdef double pi, sl_new, sr_new
    cdef double sr = 0.
    cdef double sl = 1.
    if ab * ac < a * abcd:
        for i in range(a + 1, a_max + 1):
            sr_new = sr + exp(pa - lgamma(i + 1) - lgamma(ab - i + 1) - lgamma(ac - i + 1) - lgamma(abcd - ab - ac + i + 1))
            if sr_new == sr:
                break
            sr = sr_new
        return -log(1. - max(0, exp(p0 - pa) * sr))
    else:
        for i in range(a - 1, a_min - 1, -1):
            sl_new = sl + exp(pa - lgamma(i + 1) - lgamma(ab - i + 1) - lgamma(ac - i + 1) - lgamma(abcd - ab - ac + i + 1))
            if sl_new == sl:
                break
            sl = sl_new
        return max(0, pa - p0 - log(sl))

cdef inline double mlog10Test1l(long long a, long long b, long long c, long long d) except *:
    return mlnTest2l(a, a + b, a + c, a + b + c + d) / LN10

cdef inline double mlog10Test2l(long long a, long long ab, long long ac, long long abcd) except *:
    return mlnTest2l(a, ab, ac, abcd) / LN10

# ======================== Right Tail Only ========================
cpdef double test1r(long long a, long long b, long long c, long long d) except *:
    return exp(-mlnTest2r(a, a + b, a + c, a + b + c + d))

cpdef double test2r(long long a, long long ab, long long ac, long long abcd) except *:
    return exp(-mlnTest2r(a, ab, ac, abcd))

cdef inline double mlnTest1r(long long a, long long b, long long c, long long d) except *:
    return mlnTest2r(a, a + b, a + c, a + b + c + d)

cdef inline double mlnTest2r(long long a, long long ab, long long ac, long long abcd) except *:
    if 0 > a or a > ab or a > ac or ab + ac > abcd + a:
        raise ValueError('invalid contingency table')
    if abcd > MAXN:
        raise OverflowError('the grand total of contingency table is too large')
    cdef long long a_min = max(0, ab + ac - abcd)
    cdef long long a_max = min(ab, ac)
    if a_min == a_max:
        return 0.
    cdef double p0 = lgamma(ab + 1) + lgamma(ac + 1) + lgamma(abcd - ac + 1) + lgamma(abcd - ab + 1) - lgamma(abcd + 1)
    cdef double pa = lgamma(a + 1) + lgamma(ab - a + 1) + lgamma(ac - a + 1) + lgamma(abcd - ab - ac + a + 1)
    cdef double sl_new, sr_new
    cdef double sl = 0.
    cdef double sr = 1.
    if ab * ac > a * abcd:
        for i in range(a - 1, a_min - 1, -1):
            sl_new = sl + exp(pa - lgamma(i + 1) - lgamma(ab - i + 1) - lgamma(ac - i + 1) - lgamma(abcd - ab - ac + i + 1))
            if sl_new == sl:
                break
            sl = sl_new
        return -log(1. - max(0, exp(p0 - pa) * sl))
    else:
        for i in range(a + 1, a_max + 1):
            sr_new = sr + exp(pa - lgamma(i + 1) - lgamma(ab - i + 1) - lgamma(ac - i + 1) - lgamma(abcd - ab - ac + i + 1))
            if sr_new == sr:
                break
            sr = sr_new
        return max(0, pa - p0 - log(sr))

cdef inline double mlog10Test1r(long long a, long long b, long long c, long long d) except *:
    return mlnTest2r(a, a + b, a + c, a + b + c + d) / LN10

cdef inline double mlog10Test2r(long long a, long long ab, long long ac, long long abcd) except *:
    return mlnTest2r(a, ab, ac, abcd) / LN10

# ======================== Two Tails Only ========================
cpdef double test1t(long long a, long long b, long long c, long long d) except *:
    return exp(-mlnTest2t(a, a + b, a + c, a + b + c + d))

cpdef double test2t(long long a, long long ab, long long ac, long long abcd) except *:
    return exp(-mlnTest2t(a, ab, ac, abcd))

cdef inline double mlnTest1t(long long a, long long b, long long c, long long d) except *:
    return mlnTest2t(a, a + b, a + c, a + b + c + d)

cdef inline double mlnTest2t(long long a, long long ab, long long ac, long long abcd) except *:
    if 0 > a or a > ab or a > ac or ab + ac > abcd + a:
        raise ValueError('invalid contingency table')
    if abcd > MAXN:
        raise OverflowError('the grand total of contingency table is too large')
    cdef long long a_min = max(0, ab + ac - abcd)
    cdef long long a_max = min(ab, ac)
    if a_min == a_max:
        return 0.
    cdef double p0 = lgamma(ab + 1) + lgamma(ac + 1) + lgamma(abcd - ac + 1) + lgamma(abcd - ab + 1) - lgamma(abcd + 1)
    cdef double pa = lgamma(a + 1) + lgamma(ab - a + 1) + lgamma(ac - a + 1) + lgamma(abcd - ab - ac + a + 1)
    cdef double st = 1.
    cdef double pi, st_new
    if ab * ac < a * abcd:
        for i in range(min(a - 1, llround(ab * ac / abcd)), a_min - 1, -1):
            pi = lgamma(i + 1) + lgamma(ab - i + 1) + lgamma(ac - i + 1) + lgamma(abcd - ab - ac + i + 1)
            if pi < pa:
                continue
            st_new = st + exp(pa - pi)
            if st_new == st:
                break
            st = st_new
        for i in range(a + 1, a_max + 1):
            pi = lgamma(i + 1) + lgamma(ab - i + 1) + lgamma(ac - i + 1) + lgamma(abcd - ab - ac + i + 1)
            st_new = st + exp(pa - pi)
            if st_new == st:
                break
            st = st_new
    else:
        for i in range(a - 1, a_min - 1, -1):
            pi = lgamma(i + 1) + lgamma(ab - i + 1) + lgamma(ac - i + 1) + lgamma(abcd - ab - ac + i + 1)
            st_new = st + exp(pa - pi)
            if st_new == st:
                break
            st = st_new
        for i in range(max(a + 1, llround(ab * ac / abcd)), a_max + 1):
            pi = lgamma(i + 1) + lgamma(ab - i + 1) + lgamma(ac - i + 1) + lgamma(abcd - ab - ac + i + 1)
            if pi < pa:
                continue
            st_new = st + exp(pa - pi)
            if st_new == st:
                break
            st = st_new
    return max(0, pa - p0 - log(st))

cpdef inline double mlog10Test1t(long long a, long long b, long long c, long long d) except *:
    return mlnTest2t(a, a + b, a + c, a + b + c + d) / LN10

cdef inline double mlog10Test2t(long long a, long long ab, long long ac, long long abcd) except *:
    return mlnTest2t(a, ab, ac, abcd) / LN10

cpdef double fisher_exact(long long a, long long b, long long c, long long d, str alternative: str) except *:
    """
    Perform a Fisher exact test on a 2x2 contingency table.

    :param a: row 1 col 1
    :param b: row 1 col 2
    :param c: row 2 col 1
    :param d: row 2 col 2
    :param alternative: {‘two-sided’, ‘less’, ‘greater’} (cpdefault: 'two-sided')
    :return: pvalue
    """
    if alternative == 'two-sided':
        return test1t(a, b, c, d)
    elif alternative == 'less':
        return test1l(a, b, c, d)
    elif alternative == 'greater':
        return test1r(a, b, c, d)
    else:
        raise ValueError("`alternative` should be one of {'two-sided', 'less', 'greater'}")

cpdef double odds_ratio(double a, double b, double c, double d) except *:
    """
    Calculate odds ratio of a contingency table.

    :param a: row 1 col 1
    :param b: row 1 col 2
    :param c: row 2 col 1
    :param d: row 2 col 2
    :return: odds ratio (this is prior odds ratio and not a posterior estimate.)
    """
    if a + b == 0 or c + d == 0 or a + c == 0 or b + d == 0:
        return pymath_nan

    if not (b > 0 and c > 0):
        return pymath_inf

    return (a * d) / (c * b)
