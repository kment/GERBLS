"""
setup.py for installing GERBLS

Most of the project settings/metadata are in `pyproject.toml`.
Currently, this file is included because the setuptools (v80) support for defining extention modules
inside `pyproject.toml` is still considered experimental and may be subject to change.
"""
from setuptools import Extension, setup
from Cython.Build import cythonize
import numpy as np

setup(name = "gerbls",
      ext_modules = cythonize([Extension("gerbls.core",
                                         ["src/gerbls/core.pyx",
                                          "src/gerbls/cpp/ffafunc.cpp",
                                          "src/gerbls/cpp/model.cpp",
                                          "src/gerbls/cpp/physfunc.cpp",
                                          "src/gerbls/cpp/structure.cpp"],
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
                              annotate=False,
                              compiler_directives={"embedsignature": True}
                              ),
      zip_safe = False,
      )
