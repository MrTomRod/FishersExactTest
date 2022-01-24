## Credit

This is heavily based on the work of [painyeph](https://github.com/painyeph/FishersExactTest).

I basically added numba compilation and easier installation. Nevertheless, it makes a big difference speed-wise!

## Installation

```bash
# First, python-dev / python-devel appears to be required.
# On Fedora:
sudo dnf install python3-devel
# On Debian/Ubuntu:
sudo apt-get install python3-dev

# Then, install the package using pip:
pip install git+https://github.com/MrTomRod/fast-fisher
# or
pip install fast-fisher  # from https://pypi.org/project/fast-fisher/
```

## Usage

**Recommended usage:**

```python
from fast_fisher import fast_fisher_exact, odds_ratio

a, b, c, d = 15, 32, 25, 46

pvalue = fast_fisher_exact(a, b, c, d, alternative='two-sided')
odds = odds_ratio(a, b, c, d)
```

**Alternative usages:**

`fast_fisher_exact_compatibility` has the same syntax as `scipy.stats.fisher_exact`.

```python
from math import isclose
from scipy.stats import fisher_exact
from fast_fisher import fast_fisher_exact_compatibility

table = [[15, 32], [25, 46]]

for alternative in ['two-sided', 'less', 'greater']:
    odds_s, pval_s = fisher_exact(table, alternative)
    odds_f, pval_f = fast_fisher_exact_compatibility(table, alternative)
    assert odds_s == odds_f
    assert isclose(pval_s, pval_f)  # not always true!
```

### Advanced Usage

| test type    | p-value                                                | -log( p-value )                                              | -log10( p-value )                                                  |
|--------------|--------------------------------------------------------|--------------------------------------------------------------|--------------------------------------------------------------------|
| left-tailed  | `test1l(a, b, c, d)` or `test2l(a, a+b, a+c, a+b+c+d)` | `mlnTest1l(a, b, c, d)` or `mlnTest2l(a, a+b, a+c, a+b+c+d)` | `mlog10Test1l(a, b, c, d)` or `mlog10Test2l(a, a+b, a+c, a+b+c+d)` |
| right-tailed | `test1r(a, b, c, d)` or `test2r(a, a+b, a+c, a+b+c+d)` | `mlnTest1r(a, b, c, d)` or `mlnTest2r(a, a+b, a+c, a+b+c+d)` | `mlog10Test1r(a, b, c, d)` or `mlog10Test2r(a, a+b, a+c, a+b+c+d)` |
| two-tailed   | `test1t(a, b, c, d)` or `test2t(a, a+b, a+c, a+b+c+d)` | `mlnTest1t(a, b, c, d)` or `mlnTest2t(a, a+b, a+c, a+b+c+d)` | `mlog10Test1t(a, b, c, d)` or `mlog10Test2t(a, a+b, a+c, a+b+c+d)` |
| all          | `test1(a, b, c, d)` or `test2(a, a+b, a+c, a+b+c+d)`   | `mlnTest1(a, b, c, d)` or `mlnTest2(a, a+b, a+c, a+b+c+d)`   | `mlog10Test1(a, b, c, d)` or `mlog10Test2(a, a+b, a+c, a+b+c+d)`   |

## Speed

Comparison of

- `scipy.stats.fisher_exact` ([scipy docs](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.fisher_exact.html))
- `fast_fisher.fast_fisher_python` (this library)
- `fast_fisher.fast_fisher_compiled` (this library, compiled using numba)
- `fast_fisher.fast_fisher_cpython` (this library, compiled using cython)
- `fisher.pvalue` (from [brentp/fishers_exact_test](https://github.com/brentp/fishers_exact_test))

See `benchmark.py`.

|      a |      b |      c |      d |    test type |     scipy |  f_python | f_compiled |   f_cython |     brentp |
|-------:|-------:|-------:|-------:|-------------:|----------:|----------:|-----------:|-----------:|-----------:|
|      8 |      2 |      1 |      5 |  left-tailed |    130 us |      3 us |       0 us |       0 us |       0 us |
|      8 |      2 |      1 |      5 | right-tailed |    124 us |      3 us |       0 us |       0 us |       0 us |
|      8 |      2 |      1 |      5 |   two-tailed |    895 us |      6 us |       1 us |       1 us |       0 us |
|     10 |    100 |     10 |    100 |  left-tailed |    143 us |     10 us |       1 us |       1 us |       0 us |
|     10 |    100 |     10 |    100 | right-tailed |    147 us |     12 us |       1 us |       1 us |       1 us |
|     10 |    100 |     10 |    100 |   two-tailed |    198 us |     18 us |       2 us |       2 us |       1 us |
|     10 |   1000 |  10000 | 100000 |  left-tailed |    139 us |     11 us |       1 us |       1 us |     110 us |
|    100 |   1000 |  10000 | 100000 | right-tailed |    238 us |     78 us |       6 us |       6 us |     142 us |
|    100 |   1000 |  10000 | 100000 |   two-tailed |    186 us |    137 us |      11 us |      11 us |     136 us |
|  10000 |    100 |   1000 | 100000 |  left-tailed |   1010 us |      8 us |       1 us |       1 us |    1486 us |
|  10000 |    100 |   1000 | 100000 | right-tailed |    150 us |      6 us |       1 us |       1 us |    1495 us |
|  10000 |    100 |   1000 | 100000 |   two-tailed |  63192 us |    768 us |      55 us |      58 us |    1459 us |
|  10000 |  10000 |  10000 |  10000 |  left-tailed |    967 us |    338 us |      27 us |      28 us |    2808 us |
|  10000 |  10000 |  10000 |  10000 | right-tailed |    969 us |    344 us |      27 us |      28 us |    2820 us |
|  10000 |  10000 |  10000 |  10000 |   two-tailed |    177 us |    689 us |      54 us |      61 us |    2952 us |

While numba seems be marginally faster than the cython, cython is a much better build dependency. Therefore, by default, only the cython version is
installed.

## Precision

```python
from numpy import log10, isinf
from scipy.stats import fisher_exact
from fast_fisher import fast_fisher

scipy_fisher = lambda t: fisher_exact([[t[0], t[1]], [t[2], t[3]]])[1]

print(f"{'contingency table':<30} {'fast -log10(pvalue)':<21} {'scipy -log10(pvalue)'}")
for exponent in range(0, 16):
    table = (100, 1, 10, 10 ** exponent)
    fast_mlog = fast_fisher.mlog10Test1t(*table)
    scipy_mlog = -log10(scipy_fisher(table))

    if isinf(scipy_mlog):
        scipy_mlog = 'failed to compute'

    print(f"{str(table):<30} {fast_mlog:<21} {scipy_mlog}")
```

```text
contingency table              fast -log10(pvalue)   scipy -log10(pvalue)
(100, 1, 10, 1)                0.7268124553699258    0.7268124553698625
(100, 1, 10, 10)               7.831294376070296     7.831294376070258
(100, 1, 10, 100)              46.49556750272154     46.4955675027216
(100, 1, 10, 1000)             128.93472935802276    128.93472935802373
(100, 1, 10, 10000)            226.62104816785       226.62104816785057
(100, 1, 10, 100000)           326.3812661048001     failed to compute
(100, 1, 10, 1000000)          426.35719912501844    failed to compute
(100, 1, 10, 10000000)         526.3547915160074     failed to compute
(100, 1, 10, 100000000)        626.3545507841483     failed to compute
(100, 1, 10, 1000000000)       726.3545273226812     failed to compute
(100, 1, 10, 10000000000)      826.3545146294285     failed to compute
(100, 1, 10, 100000000000)     926.354158998833      failed to compute
(100, 1, 10, 1000000000000)    1026.3583095975994    failed to compute
(100, 1, 10, 10000000000000)   1126.3427388160835    failed to compute
(100, 1, 10, 100000000000000)  1226.447616894783     failed to compute
<input>:11: RuntimeWarning: divide by zero encountered in log10
Traceback (most recent call last):
  File "/usr/lib64/python3.10/code.py", line 90, in runcode
    exec(code, self.locals)
  File "<input>", line 10, in <module>
  File "fast_fisher/fast_fisher_cython.pyx", line 236, in fast_fisher.fast_fisher_cython.mlog10Test1t
    cpdef inline double mlog10Test1t(long long a, long long b, long long c, long long d) except *:
  File "fast_fisher/fast_fisher_cython.pyx", line 237, in fast_fisher.fast_fisher_cython.mlog10Test1t
    return mlnTest2t(a, a + b, a + c, a + b + c + d) / LN10
  File "fast_fisher/fast_fisher_cython.pyx", line 195, in fast_fisher.fast_fisher_cython.mlnTest2t
    raise OverflowError('the grand total of contingency table is too large')
OverflowError: the grand total of contingency table is too large
```
