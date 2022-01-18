## Credit

This is heavily based on the work of [painyeph](https://github.com/painyeph/FishersExactTest).

I basically added numba compilation and easier installation.
Nevertheless, it makes a big difference speed-wise!

## Installation

```bash
pip install git+https://github.com/MrTomRod/FishersExactTest
```

## Functions

| test type    | p-value                                                | -log( p-value )                                              | -log10( p-value )                                                  |
|--------------|--------------------------------------------------------|--------------------------------------------------------------|--------------------------------------------------------------------|
| left-tailed  | `test1l(a, b, c, d)` or `test2l(a, a+b, a+c, a+b+c+d)` | `mlnTest1l(a, b, c, d)` or `mlnTest2l(a, a+b, a+c, a+b+c+d)` | `mlog10Test1l(a, b, c, d)` or `mlog10Test2l(a, a+b, a+c, a+b+c+d)` |
| right-tailed | `test1r(a, b, c, d)` or `test2r(a, a+b, a+c, a+b+c+d)` | `mlnTest1r(a, b, c, d)` or `mlnTest2r(a, a+b, a+c, a+b+c+d)` | `mlog10Test1r(a, b, c, d)` or `mlog10Test2r(a, a+b, a+c, a+b+c+d)` |
| two-tailed   | `test1t(a, b, c, d)` or `test2t(a, a+b, a+c, a+b+c+d)` | `mlnTest1t(a, b, c, d)` or `mlnTest2t(a, a+b, a+c, a+b+c+d)` | `mlog10Test1t(a, b, c, d)` or `mlog10Test2t(a, a+b, a+c, a+b+c+d)` |
| all          | `test1(a, b, c, d)` or `test2(a, a+b, a+c, a+b+c+d)`   | `mlnTest1(a, b, c, d)` or `mlnTest2(a, a+b, a+c, a+b+c+d)`   | `mlog10Test1(a, b, c, d)` or `mlog10Test2(a, a+b, a+c, a+b+c+d)`   |

## Speed

Comparison of `scipy.stats.fisher_exact`, `fast_fisher.fast_fisher_python` and `fast_fisher.fast_fisher_compiled`. See `benchmark.py`.

|      a |      b |      c |      d |    test type |     scipy |  f_python | f_compiled |
|-------:|-------:|-------:|-------:|-------------:|----------:|----------:|-----------:|
|      8 |      2 |      1 |      5 |  left-tailed |    130 us |      3 us |       1 us |
|      8 |      2 |      1 |      5 | right-tailed |    134 us |      3 us |       0 us |
|      8 |      2 |      1 |      5 |   two-tailed |    884 us |      6 us |       1 us |
|    100 |   1000 |  10000 | 100000 |  left-tailed |    153 us |     57 us |       5 us |
|    100 |   1000 |  10000 | 100000 | right-tailed |    218 us |     73 us |       6 us |
|    100 |   1000 |  10000 | 100000 |   two-tailed |    174 us |    130 us |      10 us |
|  10000 |    100 |   1000 | 100000 |  left-tailed |    955 us |      8 us |       1 us |
|  10000 |    100 |   1000 | 100000 | right-tailed |    135 us |      7 us |       1 us |
|  10000 |    100 |   1000 | 100000 |   two-tailed |  62821 us |    753 us |      53 us |
|  10000 |  10000 |  10000 |  10000 |  left-tailed |    900 us |    323 us |      26 us |
|  10000 |  10000 |  10000 |  10000 | right-tailed |    912 us |    329 us |      27 us |
|  10000 |  10000 |  10000 |  10000 |   two-tailed |    172 us |    679 us |      53 us |

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
Traceback (most recent call last):
  File "/usr/lib64/python3.10/code.py", line 90, in runcode
    exec(code, self.locals)
  File "<input>", line 10, in <module>
  File "/tmp/pip-req-build-usaxopn1/fast_fisher/fast_fisher_numba.py", line 244, in mlnTest2t
OverflowError: the grand total of contingency table is too large
```
