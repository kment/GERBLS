from setuptools import Extension, setup
from Cython.Build import cythonize
import numpy as np

setup(name = "GERBLS",
      ext_modules = cythonize([Extension("GERBLS", 
                                         ["GERBLS/gerbls.pyx",
                                          "GERBLS/source/ffafunc.cpp",
                                          "GERBLS/source/model.cpp",
                                          "GERBLS/source/physfunc.cpp",
                                          "GERBLS/source/structure.cpp"],
                                         include_dirs=[np.get_include()],
                                         define_macros=[("NPY_NO_DEPRECATED_API", 
                                                         "NPY_1_7_API_VERSION")],
                                         extra_compile_args = ["-O3",
                                                               "-std=c++0x",
                                                               "-march=native",
                                                               "-fassociative-math",
                                                               "-fno-math-errno",
                                                               "-ffinite-math-only",
                                                               "-fno-rounding-math",
                                                               "-fno-signed-zeros",
                                                               "-fno-trapping-math"])], 
                              annotate=False),
      zip_safe = False,
      )
