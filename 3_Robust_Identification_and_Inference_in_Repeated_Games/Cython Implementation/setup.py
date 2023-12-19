from setuptools import setup
from Cython.Build import cythonize

import Cython.Compiler.Options
Cython.Compiler.Options.annotate = True

import numpy

setup(
    ext_modules = cythonize("RIIRG_lib.pyx", annotate = True), 
    include_dirs=[numpy.get_include()]
)

# python setup.py build_ext --inplace --force