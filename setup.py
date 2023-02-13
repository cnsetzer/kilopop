from numpy.distutils.core import setup, Extension
import os
from setuptools import find_packages
import numpy.f2py as npf2py

assert "LD_LIBRARY_PATH" in os.environ.keys(), "Did you add lapack to LD_LIBRARY_PATH?"

ext = Extension(
    name="macronova2py",
    sources=[
        "bnspopkne/fortran_source/module_physics_constants.f90",
        "bnspopkne/fortran_source/hratelib.f90",
        "bnspopkne/fortran_source/macronova_Pinto_Eastman_CNS.f90",
        "bnspopkne/fortran_source/macronova2py.f90",
    ],
    extra_f90_compile_args=[
        "-cpp",
        "-g",
        "-O3",
        "-ffpe-trap=overflow,underflow,invalid",
        "-Wall",
        "-fbacktrace",
        "-fimplicit-none",
        "-fdefault-double-8",
        "-fdefault-real-8",
        "-fopenmp",
        "-ffree-line-length-512",
    ],
    libraries=["lapack", "blas"],
    library_dirs=[
        os.environ["LD_LIBRARY_PATH"],
    ],
    f2py_options=["c", "only:", "calculate_luminosity", ":", "m"],
)


setup(ext_modules=[ext])
