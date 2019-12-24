# Prerequisites

In this section, the library's hardware and software requirements are defined.

## Hardware Requirements

This library targets x86\_64 chips and, in particular, Intel chips.  The library can additionally leverage high-end NVidia GPU's.

## Software Requirements

The following packages are required:

   1.  A C++17 compliant compiler.
   2.  [CMake](https://cmake.org/) v.3.10 or greater.
   3.  [Intel Math Kernel Library](https://software.intel.com/en-us/mkl)
   4.  [Intel Performance Primitives](https://software.intel.com/en-us/ipp)
   5.  [Google Test](https://github.com/google/googletest)

The following packages are optional:

   1.  [CUDA 10.1](https://developer.nvidia.com/cuda-downloads) for a GPU-based matched filtered implementation.  I may relax this and rewrite this in OpenCL or, for my convenience, rewrite this in Data Parallel C++.
   2.  [OpenMP](https://www.openmp.org/).  This likely probably ships with your compiler.
   3.  Message Passing Interface v3, e.g., [MPICH](https://www.mpich.org/) or [OpenMPI](https://www.open-mpi.org/).
   4.  The Parallel Standard Template Library extensions.  This will be better defined in C++20.
   5.  [pybind11](https://github.com/pybind/pybind11) for generation of Python bindings.
 
# Downloading the Software

    git clone https://github.com/uofuseismo/mflib.git

# Configure the Software

After the software is downloaded, CMake must be configured as to generate a Makefile for your system.  To introduce this step I'll attempt to show off some configuration scripts.  

The simplest example assumes all prerequisites are installed in a sensible place.  In that case, the following configuration script run in the root source directory may be enough to get CMake off the ground

    #!/bin/bash
    export GTEST_ROOT=/usr/local
    export BUILD_DIR=gcc_build
    if [ -d ${BUILD_DIR} ]; then
       rm -rf ${BUILD_DIR}
    fi
    mkdir ${BUILD_DIR}
    cd ${BUILD_DIR}
    cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_FLAGS="-g -Wall -std=c++17" \
    cd ../

When available, I prefer using the Intel C++ compiler 

    #!/bin/bash
    export CXX=/opt/intel/bin/icpc
    export GTEST_ROOT=/usr/local
    export BUILD_DIR=intel_build
    if [ -d ${BUILD_DIR} ]; then
       rm -rf ${BUILD_DIR}
    fi
    mkdir ${BUILD_DIR}
    cd ${BUILD_DIR}
    cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_COMPILER=${CXX} \
    -DCMAKE_CXX_FLAGS="-g -Wall -std=c++17 -mkl=sequential -xCORE-AVX512" \
    -DMPI_CXX_INCLUDE_PATH=${MPI_ROOT_DIR}/intel64/include \
    -DMPI_CXX_LIBRARIES="/opt/intel/impi/2019.5.281/intel64/lib/release/libmpi.so;/opt/intel/compilers_and_libraries_2019.5.281/linux/mpi/intel64/libfabric//lib/libfabric.so"
    cd ../

For the University of Utah folks we have to look for things in non-standard locations.  You may need to do the following on a CHPC landing node

    #!/bin/bash
    module load gcc/8.3.0
    module load cmake/3.15.3
    module load cuda/10.1
    export CXX=/uufs/chpc.utah.edu/sys/installdir/gcc/8.3.0/bin/g++
    export MKL_ROOT=/uufs/chpc.utah.edu/sys/installdir/intel/compilers_and_libraries_2019.5.281/linux/mkl
    export IPP_ROOT=/uufs/chpc.utah.edu/sys/installdir/intel/compilers_and_libraries_2019.5.281/linux/ipp
    export GTEST_ROOT=/uufs/chpc.utah.edu/common/home/u6024448/Codes/googletest
    export CMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH}:/uufs/chpc.utah.edu/common/home/koper-group1/bbaker/templateMatchingSource/pybind11/share/cmake/pybind11
    export BUILD_DIR=gcc_build
    if [ -d ${BUILD_DIR} ]; then
       rm -rf ${BUILD_DIR}
    fi
    mkdir ${BUILD_DIR}
    cd ${BUILD_DIR}
    cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_COMPILER=${CXX} \
    -DCMAKE_CXX_FLAGS="-g -Wall -Wextra" \
    -DMFLIB_WRAP_PYTHON=TRUE

Provided this is successful then descend into the build directory and type

    make

Following a successful build you should verify the software works as advertised by typing

    make test
 
