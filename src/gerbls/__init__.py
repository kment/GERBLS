# GERBLS version
__version__ = "0.6.0"

# Compiled Cython library
from _gerbls import *

# Core GERBLS functionality
from .blsfunc import run_bls

# Optional extras
try:
    from .clean import clean_savgol
except ImportError:
    print("Warning: GERBLS has been installed without extras. Only the core BLS functionality is supported.")
    pass