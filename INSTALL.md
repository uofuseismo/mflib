# Prerequisites

## Hardware Requirements

This library targets x86\_64 chips and, in particular, Intel chips.  The library can additionally leverage high-end NVidia GPU's.

## Software Requirements

The following packages are required:

   1.  A C++17 compliant compiler with the PSTL extensions.
   2.  [CMake](https://cmake.org/) v.3.10 or greater.
   3.  [Intel Math Kernel Library](https://software.intel.com/en-us/mkl)
   4.  [Intel Performance Primitives](https://software.intel.com/en-us/ipp)
   5.  Message Passing Interface v3, e.g., [MPICH](https://www.mpich.org/) or [OpenMPI](https://www.open-mpi.org/).
   6.  [Google Test](https://github.com/google/googletest)

The following packages are optional:

   1.  [CUDA 10.1](https://developer.nvidia.com/cuda-downloads) for a GPU-based matched filtered implementation.  I may relax this and rewrite this in OpenCL.
   2.  [OpenMP](https://www.openmp.org/).  This probably ships with your compiler.
   3.  [pybind11](https://github.com/pybind/pybind11) to generate Python bindings.
 
# Downloading the Software

    git clone https://github.com/uofuseismo/mflib.git

# Configure the Software

If things are installed in a sensible place then the following configuration script may be enough to get CMake off the ground

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
    -DCMAKE_CXX_FLAGS="-g -Wall -std=c++17 -mkl=sequential -xCORE-AVX512"
    cd ../

Provided this is successful then descend into the build directory and type

    make

Following a successful build you can type

    make test
 