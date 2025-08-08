# cython: language_level = 3str
# distutils: language = c++

from gerbls.cgerbls cimport *
from libc.stdlib cimport srand
from libc.time cimport time as ctime
from libcpp cimport bool as bool_t
from libcpp.vector cimport vector
import numpy as np
cimport numpy as np

# Initialize random number generator
srand(ctime(NULL))

# These are just text inclusions
include "blsmodel.pxi"
include "struct.pxi"
include "utils.pxi"