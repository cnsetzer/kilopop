from numpy.distutils.core import setup, Extension
import os
from setuptools import find_packages
import numpy.f2py as npf2py

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
    ],
    libraries=["lapack", "blas"],
    # libraries=["scalapack", "openblas"],
    library_dirs=[
        # os.environ["SCALAPACK_LIB"],
        # os.environ["OPENBLAS_LIB"],
        "/Users/cnsetzer/software/lib/lapack-3.8.0/",
        # "/usr/lib64/",
    ],
    f2py_options=["c", "only:", "calculate_luminosity", ":", "m"],
)


if __name__ == "__main__":

    setup(
        name="bnspopkne",
        version="0.0.1",
        description="Kilonova population model used in C. N. Setzer et al. 2022.",
        url="http://github.com/cnsetzer/Setzer2022_BNSpopkNe",
        author="Christian Setzer",
        author_email="christian.setzer@fysik.su.se",
        license="MIT",
        classifiers=[
            # How mature is this project? Common values are
            #   3 - Alpha
            #   4 - Beta
            #   5 - Production/Stable
            "Development Status :: 3 - Alpha",
            # Indicate who your project is intended for
            "Intended Audience :: Astronomers",
            # Pick your license as you wish (should match "license" above)
            "License :: OSI Approved :: MIT License",
            # Specify the Python versions you support here. In particular, ensure
            # that you indicate whether you support Python 2, Python 3 or both.
            "Programming Language :: Python :: 3.6+",
        ],
        packages=find_packages(),
        ext_modules=[ext],
    )