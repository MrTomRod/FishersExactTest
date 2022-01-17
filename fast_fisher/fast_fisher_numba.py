#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: eph, Thomas Roder

from numba.pycc import CC

cc = CC('fast_fisher_compiled')

from math import log, exp

from math import lgamma

from numba import njit, jit_module


def _maxn():
    l, n, h = 1, 2, float('inf')
    while l < n:
        if abs(lgamma(n + 1) - lgamma(n) - log(n)) >= 1:
            h = n
        else:
            l = n
        n = (l + min(h, l * 3)) // 2
    return n


LN10 = log(10)
NINF = float('-inf')
MAXN = _maxn()


# ======================== Full Test ========================


@cc.export('test1', '(i8, i8, i8, i8)')
def test1(a, b, c, d):
    result = mlnTest2(a, a + b, a + c, a + b + c + d)
    return exp(-result[0]), exp(-result[1]), exp(-result[2])


@cc.export('test2', '(i8, i8, i8, i8)')
def test2(a, ab, ac, abcd):
    result = mlnTest2(a, ab, ac, abcd)
    return exp(-result[0]), exp(-result[1]), exp(-result[2])


@cc.export('mlnTest1', '(i8, i8, i8, i8)')
def mlnTest1(a, b, c, d):
    return mlnTest2(a, a + b, a + c, a + b + c + d)


@cc.export('mlnTest2', '(i8, i8, i8, i8)')
def mlnTest2(a, ab, ac, abcd):
    if 0 > a or a > ab or a > ac or ab + ac > abcd + a:
        raise ValueError('invalid contingency table')
    if abcd > MAXN:
        raise OverflowError('the grand total of contingency table is too large')
    a_min = max(0, ab + ac - abcd)
    a_max = min(ab, ac)
    if a_min == a_max:
        return 0., 0., 0.
    p0 = lgamma(ab + 1) + lgamma(ac + 1) + lgamma(abcd - ac + 1) + lgamma(abcd - ab + 1) - lgamma(abcd + 1)
    pa = lgamma(a + 1) + lgamma(ab - a + 1) + lgamma(ac - a + 1) + lgamma(abcd - ab - ac + a + 1)
    sl = sr = 0.
    if ab * ac < a * abcd:
        for i in range(min(a - 1, int(round(ab * ac / abcd))), a_min - 1, -1):
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
        for i in range(max(a + 1, int(round(ab * ac / abcd))), a_max + 1):
            pi = lgamma(i + 1) + lgamma(ab - i + 1) + lgamma(ac - i + 1) + lgamma(abcd - ab - ac + i + 1)
            if pi < pa:
                continue
            sr_new = sr + exp(pa - pi)
            if sr_new == sr:
                break
            sr = sr_new
        return max(0, pa - p0 - log(sl + 1.)), -log(1. - max(0, exp(p0 - pa) * sl)), max(0, pa - p0 - log(sl + 1. + sr))


@cc.export('mlog10Test1', '(i8, i8, i8, i8)')
def mlog10Test1(a, b, c, d):
    result = mlnTest2(a, a + b, a + c, a + b + c + d)
    return result[0] / LN10, result[1] / LN10, result[2] / LN10


@cc.export('mlog10Test2', '(i8, i8, i8, i8)')
def mlog10Test2(a, ab, ac, abcd):
    result = mlnTest2(a, ab, ac, abcd)
    return result[0] / LN10, result[1] / LN10, result[2] / LN10


# ======================== Left Tail Only ========================
@cc.export('test1l', 'f8(i8, i8, i8, i8)')
def test1l(a, b, c, d):
    return exp(-mlnTest2l(a, a + b, a + c, a + b + c + d))


@cc.export('test2l', 'f8(i8, i8, i8, i8)')
def test2l(a, ab, ac, abcd):
    return exp(-mlnTest2l(a, ab, ac, abcd))


@cc.export('mlnTest1l', '(i8, i8, i8, i8)')
def mlnTest1l(a, b, c, d):
    return mlnTest2l(a, a + b, a + c, a + b + c + d)


@cc.export('mlnTest2l', 'f8(i8, i8, i8, i8)')
def mlnTest2l(a, ab, ac, abcd):
    if 0 > a or a > ab or a > ac or ab + ac > abcd + a:
        raise ValueError('invalid contingency table')
    if abcd > MAXN:
        raise OverflowError('the grand total of contingency table is too large')
    a_min = max(0, ab + ac - abcd)
    a_max = min(ab, ac)
    if a_min == a_max:
        return 0.
    p0 = lgamma(ab + 1) + lgamma(ac + 1) + lgamma(abcd - ac + 1) + lgamma(abcd - ab + 1) - lgamma(abcd + 1)
    pa = lgamma(a + 1) + lgamma(ab - a + 1) + lgamma(ac - a + 1) + lgamma(abcd - ab - ac + a + 1)
    if ab * ac < a * abcd:
        sr = 0.
        for i in range(a + 1, a_max + 1):
            sr_new = sr + exp(pa - lgamma(i + 1) - lgamma(ab - i + 1) - lgamma(ac - i + 1) - lgamma(abcd - ab - ac + i + 1))
            if sr_new == sr:
                break
            sr = sr_new
        return -log(1. - max(0, exp(p0 - pa) * sr))
    else:
        sl = 1.
        for i in range(a - 1, a_min - 1, -1):
            sl_new = sl + exp(pa - lgamma(i + 1) - lgamma(ab - i + 1) - lgamma(ac - i + 1) - lgamma(abcd - ab - ac + i + 1))
            if sl_new == sl:
                break
            sl = sl_new
        return max(0, pa - p0 - log(sl))


@cc.export('mlog10Test1l', 'f8(i8, i8, i8, i8)')
def mlog10Test1l(a, b, c, d):
    return mlnTest2l(a, a + b, a + c, a + b + c + d) / LN10


@cc.export('mlog10Test2l', 'f8(i8, i8, i8, i8)')
def mlog10Test2l(a, ab, ac, abcd):
    return mlnTest2l(a, ab, ac, abcd) / LN10


# ======================== Right Tail Only ========================
@cc.export('test1r', 'f8(i8, i8, i8, i8)')
def test1r(a, b, c, d):
    return exp(-mlnTest2r(a, a + b, a + c, a + b + c + d))


@cc.export('test2r', 'f8(i8, i8, i8, i8)')
def test2r(a, ab, ac, abcd):
    return exp(-mlnTest2r(a, ab, ac, abcd))


@cc.export('mlnTest1r', 'f8(i8, i8, i8, i8)')
def mlnTest1r(a, b, c, d):
    return mlnTest2r(a, a + b, a + c, a + b + c + d)


@cc.export('mlnTest2r', 'f8(i8, i8, i8, i8)')
def mlnTest2r(a, ab, ac, abcd):
    if 0 > a or a > ab or a > ac or ab + ac > abcd + a:
        raise ValueError('invalid contingency table')
    if abcd > MAXN:
        raise OverflowError('the grand total of contingency table is too large')
    a_min = max(0, ab + ac - abcd)
    a_max = min(ab, ac)
    if a_min == a_max:
        return 0.
    p0 = lgamma(ab + 1) + lgamma(ac + 1) + lgamma(abcd - ac + 1) + lgamma(abcd - ab + 1) - lgamma(abcd + 1)
    pa = lgamma(a + 1) + lgamma(ab - a + 1) + lgamma(ac - a + 1) + lgamma(abcd - ab - ac + a + 1)
    if ab * ac > a * abcd:
        sl = 0.
        for i in range(a - 1, a_min - 1, -1):
            sl_new = sl + exp(pa - lgamma(i + 1) - lgamma(ab - i + 1) - lgamma(ac - i + 1) - lgamma(abcd - ab - ac + i + 1))
            if sl_new == sl:
                break
            sl = sl_new
        return -log(1. - max(0, exp(p0 - pa) * sl))
    else:
        sr = 1.
        for i in range(a + 1, a_max + 1):
            sr_new = sr + exp(pa - lgamma(i + 1) - lgamma(ab - i + 1) - lgamma(ac - i + 1) - lgamma(abcd - ab - ac + i + 1))
            if sr_new == sr:
                break
            sr = sr_new
        return max(0, pa - p0 - log(sr))


@cc.export('mlog10Test1r', 'f8(i8, i8, i8, i8)')
def mlog10Test1r(a, b, c, d):
    return mlnTest2r(a, a + b, a + c, a + b + c + d) / LN10


@cc.export('mlog10Test2r', 'f8(i8, i8, i8, i8)')
def mlog10Test2r(a, ab, ac, abcd):
    return mlnTest2r(a, ab, ac, abcd) / LN10


# ======================== Two Tails Only ========================
@cc.export('test1t', 'f8(i8, i8, i8, i8)')
def test1t(a, b, c, d):
    return exp(-mlnTest2t(a, a + b, a + c, a + b + c + d))


@cc.export('test2t', 'f8(i8, i8, i8, i8)')
def test2t(a, ab, ac, abcd):
    return exp(-mlnTest2t(a, ab, ac, abcd))


@cc.export('mlnTest1t', 'f8(i8, i8, i8, i8)')
def mlnTest1t(a, b, c, d):
    return mlnTest2t(a, a + b, a + c, a + b + c + d)


@cc.export('mlnTest2t', 'f8(i8, i8, i8, i8)')
def mlnTest2t(a, ab, ac, abcd):
    if 0 > a or a > ab or a > ac or ab + ac > abcd + a:
        raise ValueError('invalid contingency table')
    if abcd > MAXN:
        raise OverflowError('the grand total of contingency table is too large')
    a_min = max(0, ab + ac - abcd)
    a_max = min(ab, ac)
    if a_min == a_max:
        return 0.
    p0 = lgamma(ab + 1) + lgamma(ac + 1) + lgamma(abcd - ac + 1) + lgamma(abcd - ab + 1) - lgamma(abcd + 1)
    pa = lgamma(a + 1) + lgamma(ab - a + 1) + lgamma(ac - a + 1) + lgamma(abcd - ab - ac + a + 1)
    st = 1.
    if ab * ac < a * abcd:
        for i in range(min(a - 1, int(round(ab * ac / abcd))), a_min - 1, -1):
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
        for i in range(max(a + 1, int(round(ab * ac / abcd))), a_max + 1):
            pi = lgamma(i + 1) + lgamma(ab - i + 1) + lgamma(ac - i + 1) + lgamma(abcd - ab - ac + i + 1)
            if pi < pa:
                continue
            st_new = st + exp(pa - pi)
            if st_new == st:
                break
            st = st_new
    return max(0, pa - p0 - log(st))


@cc.export('mlog10Test1t', 'f8(i8, i8, i8, i8)')
def mlog10Test1t(a, b, c, d):
    return mlnTest2t(a, a + b, a + c, a + b + c + d) / LN10


@cc.export('mlog10Test2t', 'f8(i8, i8, i8, i8)')
def mlog10Test2t(a, ab, ac, abcd):
    return mlnTest2t(a, ab, ac, abcd) / LN10


jit_module(
    nopython=True,
    cache=True,
    nogil=True
)

if __name__ == '__main__':
    cc.compile()
