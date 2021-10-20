from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy 
from Cython.Distutils import build_ext

import subprocess as sbp
import os.path as osp

import os

# Edit to your "WORK" directory
WORK = "/mwork0/sugiymnn/WORK"

COSMO = WORK + "/cosmo"

GSL_INCLUDE = COSMO + "/gsl/include"
GSL_LIB = COSMO + "/gsl/lib"

FFTW_INCLUDE = COSMO + "/fftw3/include"
FFTW_LIB     = COSMO + "/fftw3/lib"

hitomipy_ext = Extension("hitomipy", sources=["hitomi.pyx"], 
                    extra_compile_args=['-std=c++11'],
                    include_dirs=[".", numpy.get_include(), GSL_INCLUDE, FFTW_INCLUDE],
                    library_dirs=[".", GSL_LIB, FFTW_LIB],
                    libraries=['m', 'fftw3', 'gsl', 'gslcblas'],language="c++")

setup(name="hitomipy", cmdclass = {'build_ext': build_ext}, ext_modules=[hitomipy_ext])

