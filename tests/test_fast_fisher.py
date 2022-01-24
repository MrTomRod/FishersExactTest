from .utils import *


class TestFisher(TestCase):
    def _test(self, function: Callable, name: str, table, scipy_l, scipy_r, scipy_t):
        l, r, t = function(*table)

        for pval in (l, r, t):
            self.assertLessEqual(pval, 1, f'Error in {name}: {pval=} must be smaller than 1! {table=}')

        self.assertAlmostEqual(l, scipy_l, msg=f'Error in {name}: {l=} != {scipy_l}; {table=}')
        self.assertAlmostEqual(r, scipy_r, msg=f'Error in {name}: {r=} != {scipy_r}; {table=}')
        self.assertAlmostEqual(t, scipy_t, msg=f'Error in {name}: {t=} != {scipy_t}; {table=}')

    def test_with_scipy(self, samples=100):
        """
        Usually, fast-fisher and scipy agree with these randomly generated tables
        """
        for _ in range(samples):
            a, b, c, d = [int(10 ** (random() * 4) - 1) for _ in range(4)]

            scipy_l = scipy_fisher_exact(a, b, c, d, 'less')
            scipy_r = scipy_fisher_exact(a, b, c, d, 'greater')
            scipy_t = scipy_fisher_exact(a, b, c, d, 'two-sided')

            # self._test(fast_fisher_python.test1, 'fast_fisher_python', (a, b, c, d), scipy_l, scipy_r, scipy_t)
            # self._test(fast_fisher_numba.test1, 'fast_fisher_numba', (a, b, c, d), scipy_l, scipy_r, scipy_t)
            # self._test(fast_fisher_compiled.test1, 'fast_fisher_compiled', (a, b, c, d), scipy_l, scipy_r, scipy_t)
            self._test(fast_fisher_cython.test1, 'fast_fisher_compiled', (a, b, c, d), scipy_l, scipy_r, scipy_t)

    def _test_with_scipy_n(self, samples, range_max, alternative):
        n_fails = 0
        n_sig = 0
        not_identical_min_pval = 1
        for _ in range(samples):
            table = [int(random() * (range_max + 1)) for _ in range(4)]

            scipy_p = scipy_fisher_exact(*table, alternative)
            fast_p = fast_fisher_exact(*table, alternative)

            if fast_p < 0.05:
                n_sig += 1

            if not isclose(fast_p, scipy_p):
                # print(f'{table=} {fast_p < scipy_p=} {min(fast_p, scipy_p)=}')
                n_fails += 1
                not_identical_min_pval = min(not_identical_min_pval, fast_p, scipy_p)

        summary = {
            'fail_ratio': n_fails / samples,
            'significant_ratio': n_sig / samples,
            'not_identical_min_pval': not_identical_min_pval
        }
        return summary

    def test_with_scipy_n(self, samples=1000, alternative='two-sided'):
        """
        The bigger the numbers in the contingency table, the lower the disagreement between
        scipy and fast-fisher.

        range_max=10:  {'fail_ratio': 0.013, 'significant_ratio': 0.163, 'not_identical_min_pval': 0.0054179566563467355}
        range_max=20:  {'fail_ratio': 0.006, 'significant_ratio': 0.353, 'not_identical_min_pval': 4.548517044597358e-06}
        range_max=40:  {'fail_ratio': 0.005, 'significant_ratio': 0.507, 'not_identical_min_pval': 4.496376367645278e-08}
        range_max=80:  {'fail_ratio': 0.003, 'significant_ratio': 0.632, 'not_identical_min_pval': 0.4550352573371027}
        range_max=160: {'fail_ratio': 0.001, 'significant_ratio': 0.749, 'not_identical_min_pval': 1}
        """
        all_fails = {}
        for range_max in [10, 20, 40, 80, 160]:
            fails = self._test_with_scipy_n(samples, range_max, alternative)
            all_fails[range_max] = fails
            print(f'{range_max=}: {fails}')
        print(all_fails)

    def test_exception(self):
        """
        Some of these tables are examples where scipy and fast-fisher disagree.

        Notably, the compiled version of fast-fisher and the python version give somewhat differently results
        """
        for table in [
            (2, 3, 0, 2),
            (10, 13, 5, 18),
            (14, 3, 7, 18),
            (3, 14, 18, 7),
            (8, 12, 7, 13),
            (13, 12, 7, 8),
            (4, 0, 3, 7),
            (7, 3, 0, 4)
        ]:
            a, b, c, d = table
            scipy_t = scipy_fisher_exact(a, b, c, d, 'two-sided')
            # fast_t = fast_fisher_python.test1t(a, b, c, d)
            # fast_t = fast_fisher_compiled.test1t(a, b, c, d)
            fast_t = fast_fisher_cython.test1t(a, b, c, d)

            # For this contingency table, scipy and this implementation do not agree!
            try:
                self.assertAlmostEqual(fast_t, scipy_t, msg=f'{table=}')
                print(f'Agreement with {table=}')
            except AssertionError as e:
                print(f'Disagreement with {table=}', e)

    def test_swap(self, alternative='two-sided'):
        """
        Unlike in scipy's implementation, swapping can produce minor differences (~ 0.6 % in this test)
        """
        all_combinations = list(combinations(range(20), 4))
        not_identical = 0
        not_identical_min_pval = 1
        for table_orig in all_combinations:
            table_swap = fisher_swap(*table_orig)
            orig_p = fast_fisher_exact(*table_orig, alternative)
            swap_p = fast_fisher_exact(*table_swap, alternative)

            if not isclose(orig_p, swap_p):
                not_identical += 1
                not_identical_min_pval = min(not_identical_min_pval, orig_p, swap_p)

        print(f'{not_identical=} out of {len(all_combinations)}; {not_identical_min_pval=}')

    def test_ranking_preserved_naive(self, alternative='two-sided'):
        """
        conclusion: only minor disagreements at top ranks, which are what matters
        """
        res = {}
        for table_orig in combinations(range(20), 4):
            table_swap = fisher_swap(*table_orig)
            scipy_p = scipy_fisher_exact(*table_swap, alternative)
            fast_p = fast_fisher_exact(*table_swap, alternative)
            res[','.join(str(p) for p in table_swap)] = (scipy_p, fast_p)

        disagreement_df, first_disagreement_rank = get_disagreement_df(res)
        print(disagreement_df)

    def test_ranking_preserved_constant_row_sums(self, alternative='two-sided'):
        """
        conclusion: only minor disagreements at top ranks, which are what matters
        """
        for row_sum_1, row_sum_2 in [
            (20, 20),
            (15, 15),
            (10, 20),
            (5, 15),
        ]:
            rows_1 = [(i, row_sum_1 - i) for i in range(row_sum_1)]
            rows_2 = [(i, row_sum_2 - i) for i in range(row_sum_2)]
            tables = [(r1_1, r1_2, r2_1, r2_2) for (r1_1, r1_2), (r2_1, r2_2) in product(rows_1, rows_2)]
            pvals = {','.join(str(v) for v in table): (scipy_fisher_exact(*table, alternative), fast_fisher_exact(*table, alternative))
                     for table in tables}

            disagreement_df, first_disagreement_rank = get_disagreement_df(pvals)

            print('row sums: ', row_sum_1, row_sum_2)
            # print(disagreement_df)
            print('rank_scipy', first_disagreement_rank.rank_scipy, 'rank_fast', first_disagreement_rank.rank_fast)

    def test_compatibility_function(self):
        """
        Test fast_fisher_exact_compatibility, which takes the same input as scipy.stats.fisher_exact
        """

        n_not_equal = 0
        n_range = 1000
        for i in range(n_range):
            a, b, c, d = (randint(0, 50) for _ in range(4))
            or_orig, pval_orig = fisher_exact([[a, b], [c, d]])
            or_calc, pval_calc = fast_fisher_exact_compatibility([[a, b], [c, d]])

            if not isclose(pval_orig, pval_calc):
                n_not_equal += 1

            self.assertTrue(is_equivalent(or_orig, or_calc), msg=f'{or_orig=} != {or_calc=}; {[[a, b], [c, d]]}')

        print(f'{n_not_equal=} out of {n_range}')
