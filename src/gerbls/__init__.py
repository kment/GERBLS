# GERBLS version
__version__ = "0.6.3"

# Compiled Cython library
from _gerbls import *

# Core GERBLS functionality
from .blsfunc import run_bls

# Optional extras
try:
    from .clean import clean_savgol
except ImportError:
    pass