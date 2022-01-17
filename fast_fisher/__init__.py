from . import fast_fisher_python

try:
    from . import fast_fisher_compiled
    fast_fisher = fast_fisher_compiled
except ImportError as e:
    from logging import warning
    warning(f'Failed to import fast_fisher_compiled: {str(e)}')
    fast_fisher = fast_fisher_python

try:
    from . import fast_fisher_numba
except ImportError as e:
    from logging import warning
    warning(f'Failed to import fast_fisher_numba: {str(e)}')
    pass
