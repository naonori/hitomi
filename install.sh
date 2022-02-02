#!/bin/sh

#######
# Specify the directory where you want to work and the home directory where your ".bashrc" file for "export PATH" is located.
# You can put this "hitomi" directory anywhere you want.
#######
WORK=/mwork0/sugiymnn/WORK
HOME=/home/sugiymnn

#######
# Create WORK/cosmo and WORK/cosmo/temp.
# The "TEMP" directory is temporarily needed for the installation.
# All the necessary environment will be installed in WORK/cosmo.
#######
CURRENT_DIR=$(cd $(dirname $0); pwd)
ENV=$CURRENT_DIR/env
SRC=$CURRENT_DIR/src
STRUCT=$CURRENT_DIR/struct
COSMO=$WORK/cosmo
TEMP=$COSMO/temp
cd $WORK
mkdir -p $COSMO
mkdir -p $TEMP

#######
# Specifies the name of the file or directory to be used. 
# When you update the version of these libraries, modify here.
#######
name_Anaconda3=Anaconda3-2021.05-Linux-x86_64
name_fftw3=fftw-3.3.9
name_gsl=gsl-2.7
name_cuba=Cuba-4.2.1
name_cmake=cmake-3.21.1
name_lapack=lapack-3.10.0
name_multinest=MultiNest-master
name_pymultinest=PyMultiNest-2.10
name_astropy=astropy-4.3
name_class=class_public-3.0.1
name_montepython=montepython_public-3.4
name_fftlog=FFTLog-master

#######
# Install python using Anaconda3 
# Answer "yes" to all the questions.
#######
#cp -rf $ENV/$name_Anaconda3\.sh $TEMP/$name_Anaconda3\.sh
#cd $TEMP
#chmod u+x $name_Anaconda3\.sh
#Anaconda3=$COSMO/anaconda3
#./$name_Anaconda3\.sh -p $Anaconda3
#rm -rf $TEMP/$name_Anaconda3\.sh

#######
# Install fftw3
# If you want to use MPI, add "--enable-mpi" to ". /configure --prefix=$FFTW3".
# Add the "-fPIC" option to "CFLAGS" to use fftw3 in python.
#######
#cp -rf $ENV/$name_fftw3\.tar.gz $TEMP/$name_fftw3\.tar.gz
#cd $TEMP
#tar -zxvf $name_fftw3\.tar.gz
#cd $name_fftw3
#FFTW3=$COSMO/fftw3
#./configure --prefix=$FFTW3 CFLAGS="-g -O3 -fPIC" 
#make
#make install
#rm -rf $TEMP/$name_fftw3*

#######
# Install GSL
#######
#cp -rf $ENV/$name_gsl\.tar.gz $TEMP/$name_gsl\.tar.gz
#cd $TEMP
#tar -zxvf $name_gsl\.tar.gz
#cd $name_gsl
#GSL=$COSMO/gsl
#./configure --prefix=$GSL
#make
#make install
#echo "export LD_LIBRARY_PATH="$GSL/lib:\$LD_LIBRARY_PATH"" >> $HOME/.bashrc
#rm -rf $TEMP/$name_gsl*

#######
# Install Cuba.
# In order to install "pycuba" later, the "libcuba.so" file is required. 
# However, only "libcuba.a" is generated in the original Cuba library, so you have to create "libcuba.so" from "libcuba.a".
# To do so, add the "-fPIC" option to "CFLAGS" in "makefile", and install Cuba. 
# Then, use "ar" command to break "libcuba.a" into *.o files, and use "ld" command to bundle those *.o files into "libcuba.so".
# Finally, include the "cuba/lib" containing "libcuba.so" directory in your LD_LIBRARY_PATH.
#######
#cp -rf $ENV/$name_cuba\.tar.gz $TEMP/$name_cuba\.tar.gz
#cd $TEMP
#tar -zxvf $name_cuba\.tar.gz
#cd $name_cuba
#CUBA=$COSMO/cuba
#./configure --prefix=$CUBA
#sed 's/CFLAGS = -O3 -fomit-frame-pointer/CFLAGS = -O3 -fPIC -fomit-frame-pointer/g' --in-place makefile
#make
#make install
#cd $CUBA/lib
#ar xv libcuba.a
#ld -shared *.o -o libcuba.so
#rm -rf *.o
#echo "export LD_LIBRARY_PATH="$CUBA/lib:\$LD_LIBRARY_PATH"" >> $HOME/.bashrc
#rm -rf $TEMP/$name_cuba*

#######
# Install cmake.
# To install pymultinest, cmake is required. 
# If you already have cmake installed in your environment, skip this section.
#######
#cp -rf $ENV/$name_cmake\.tar.gz $TEMP/$name_cmake\.tar.gz
#cd $TEMP
#tar -zxvf $name_cmake\.tar.gz
#cd $name_cmake 
#CMAKE=$COSMO/cmake
#./bootstrap --prefix=$CMAKE
#make
#make install
#echo "export PATH="$CMAKE/bin:\$PATH"" >> $HOME/.bashrc
#rm -rf $TEMP/$name_cmake*

#######
# Install lapack.
# To install pymultinest, lapack is required. 
# If you already have lapack installed in your environment, skip this section.
# In order to generate *.so files, the "-fPIC" option needs to be added to "CMakeLists.txt file".
#######
#cp -rf $ENV/$name_lapack\.tar.gz $TEMP/$name_lapack\.tar.gz
#cd $TEMP
#tar -zxvf $name_lapack\.tar.gz
#cd $name_lapack
#LAPACK=$COSMO/lapack
#sed 's/set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DWeirdNEC -DLAPACK_ILP64 -DHAVE_LAPACK_CONFIG_H")/set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DWeirdNEC -DLAPACK_ILP64 -DHAVE_LAPACK_CONFIG_H -fPIC")/g' --in-place CMakeLists.txt
#sed 's/set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-integer-8")/set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-integer-8 -fPIC")/g' --in-place CMakeLists.txt
#sed 's/set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -recursive"/set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -recursive -fPIC"/g' --in-place CMakeLists.txt
#sed 's/set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -frecursive"/set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -frecursive -fPIC"/g' --in-place CMakeLists.txt
#sed 's/set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Mrecursive"/set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Mrecursive -fPIC"/g' --in-place CMakeLists.txt
#sed 's/set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fp-model strict")/set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fp-model strict -fPIC")/g' --in-place CMakeLists.txt
#sed 's/set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qnosave -qstrict=none")/set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qnosave -qstrict=none -fPIC")/g' --in-place CMakeLists.txt
#mkdir -p build
#cd build
#cmake -DCMAKE_INSTALL_LIBDIR=$LAPACK ..
#cmake --build . -j --target install
#echo "export LD_LIBRARY_PATH="$LAPACK:\$LD_LIBRARY_PATH"" >> $HOME/.bashrc
#rm -rf $TEMP/$name_lapack*

#######
# Install MultiNest.
#######
#cp -rf $ENV/$name_multinest\.zip $TEMP/$name_multinest\.zip
#cd $TEMP
#unzip $name_multinest\.zip
#cd $name_multinest
#cd MultiNest_v3.12_CMake/multinest
#mkdir build
#cd build
#MULTINEST=$COSMO/multinest
#cmake -DCMAKE_INSTALL_PREFIX=$MULTINEST ..
#make
#make install
#echo "export LD_LIBRARY_PATH="$MULTINEST/lib:\$LD_LIBRARY_PATH"" >> $HOME/.bashrc
#rm -rf $TEMP/$name_multinest*

#######
# Install PyMultiNest.
# Pycuba is also installed at the same time.
# Since my environment does not have access to an external network, 
# I have commented out the part in setup.py that requires connection to an external network.
# The commented out part is for running test calculations, so it does not affect the operation of the program.
#######
#cp -rf $ENV/$name_pymultinest\.zip $TEMP/$name_pymultinest\.zip
#cd $TEMP
#unzip $name_pymultinest\.zip
#cp -rf $name_pymultinest $COSMO/
#cd $COSMO
#mv $name_pymultinest pymultinest
#cd pymultinest
#sed 's/setup_requires/#setup_requires/g' --in-place setup.py
#python setup.py install
#rm -rf $TEMP/$name_pymultinest*

#######
# Install astropy.
#######
#cp -rf $ENV/$name_astropy\.zip $TEMP/$name_astropy\.zip
#cd $TEMP
#unzip $name_astropy\.zip
#cp -rf $name_astropy $COSMO/
#cd $COSMO
#mv $name_astropy astropy
#cd astropy 
#pip install astropy --no-deps
#rm -rf $TEMP/$name_astropy*

#######
# Install class
#######
#cp -rf $ENV/$name_class\.zip $TEMP/$name_class\.zip
#cd $TEMP
#unzip $name_class\.zip
#cp -rf $name_class $COSMO/
#cd $COSMO
#mv $name_class class
#cd class
#make
#cd python
#python setup.py build
#python setup.py install
#rm -rf $TEMP/$name_class*

#######
# Install montepython
# To avoid conflicts when creating your "likelihoods" directory,
# rename "montepython/montepython/likelihoods" to "montepython/montepython/likelihoods_base".
# By adding a path to the "montepython/montepython" directory, you can use "Montepython.py" anywhere.
# In montepython/run.py, add "cosmo.set(data.cosmo_arguments)" and "cosmo.compute()". 
# This is necessary in order to use the class computed for a given fiducial cosmology while varying only nuisance parameters.
#######
#cp -rf $ENV/$name_montepython\.zip $TEMP/$name_montepython\.zip
#cd $TEMP
#unzip $name_montepython\.zip
#cp -rf $name_montepython $COSMO/
#cd $COSMO
#mv $name_montepython montepython 
#cd montepython/montepython
#mv likelihoods likelihoods_base
#sed 's/sampler.run(cosmo, data, command_line)/cosmo.set(data.cosmo_arguments)\n    cosmo.compute()\n    sampler.run(cosmo, data, command_line)/g' --in-place run.py
#echo "export PYTHONPATH="$COSMO/montepython/montepython:\$PYTHONPATH"" >> $HOME/.bashrc
#rm -rf $TEMP/$name_montepython*

#######
# Install fftlog.
# Rewrite the source file of "fftlog" that you downloaded to wrap it in our python code.
#######
#cp -rf $ENV/$name_fftlog\.zip $TEMP/$name_fftlog\.zip
#cd $TEMP
#unzip $name_fftlog\.zip
#cp -rf $name_fftlog $COSMO/
#cd $COSMO
#mv $name_fftlog fftlog
#cd fftlog
#cp -rf src/fftlog.c pyfftlog.hpp
#sed 's/math.h/cmath/g' --in-place pyfftlog.hpp
#sed 's/#include <stdlib.h>/\/\/#include <stdlib.h>/g' --in-place pyfftlog.hpp
#sed 's/#include "fftlog.h"/\/\/#include "fftlog.h"/g' --in-place pyfftlog.hpp
#sed 's/I/_I_/g' --in-place pyfftlog.hpp
#sed 's/M_P_I_/M_PI/g' --in-place pyfftlog.hpp
#sed 's/EST_I_MATE/ESTIMATE/g' --in-place pyfftlog.hpp
#sed '1a\std::complex<double> _I_(0.0,1.0);' --in-place pyfftlog.hpp
#sed '1a\#include <complex>' --in-place pyfftlog.hpp
#sed 's/double complex/std::complex<double>/g' --in-place pyfftlog.hpp
#sed 's/malloc(sizeof(complex double)\*N)/new std::complex<double>\[N\]/g' --in-place pyfftlog.hpp
#sed 's/malloc (sizeof(complex double)\*N)/new std::complex<double>\[N\]/g' --in-place pyfftlog.hpp
#sed 's/creal(z)/z.real()/g' --in-place pyfftlog.hpp
#sed 's/creal(w)/w.real()/g' --in-place pyfftlog.hpp
#sed 's/cimag(w)/w.imag()/g' --in-place pyfftlog.hpp
#sed 's/creal(u\[N\/2\])/u\[N\/2\].real()/g' --in-place pyfftlog.hpp
#sed 's/creal(pow(2\*M\_PI\*r\[i\]\, -1.5) \* b\[i\])/(pow(2\*M\_PI\*r\[i\]\, -1.5) \* b\[i\]).real()/g' --in-place pyfftlog.hpp
#sed 's/free(ulocal)/delete [] ulocal/g' --in-place pyfftlog.hpp
#sed 's/free(a)/delete [] a/g' --in-place pyfftlog.hpp
#sed 's/free(b)/delete [] b/g' --in-place pyfftlog.hpp
#sed 's/cpow/pow/g' --in-place pyfftlog.hpp
#sed 's/cexp/exp/g' --in-place pyfftlog.hpp
#sed 's/clog/log/g' --in-place pyfftlog.hpp

######################
######################
######################
# From here, you can install the code that we developed.
######################
######################
######################

#######
# Copy the necessary directory structure under $WORK. 
#######
#cp -rf $STRUCT/* $WORK/

#######
# Install "hitomi_measurement".
#######
#cp -rf $SRC/hitomi_measurement $COSMO/hitomi_measurement

#######
# Install "hitomi_theory".
# Copy the "pyfftlog.hpp" file generated from the fftlog source file into "hitomi_theory/cpp".
#######
#cp -rf $SRC/hitomi_theory $COSMO/hitomi_theory
#cd $COSMO/hitomi_theory
#cp -rf $COSMO/fftlog/pyfftlog.hpp $COSMO/hitomi_theory/cpp/
#make
#echo "export PYTHONPATH="$COSMO/hitomi_theory:\$PYTHONPATH"" >> $HOME/.bashrc

rm -rf $TEMP

