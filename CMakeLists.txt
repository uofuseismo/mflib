cmake_minimum_required(VERSION 3.10)
project(mflib VERSION 0.1.2 LANGUAGES CXX)
enable_testing()

include(CheckCXXCompilerFlag)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
include(FindOpenMP QUIET)
if (OPENMP_FOUND)
   message("OpenMP found")
   string(APPEND CMAKE_CXX_FLAGS " ${OpenMP_CXX_FLAGS}")
endif()

find_package(GTest REQUIRED)
find_package(MPI)
if (MPI_FOUND)
   message("MPI was found")
endif()
find_package(CUDA
             COMPONENTS CUDA_INCLUDE_DIRS CUDA_CUFFT_LIBRARIES)
if (${CUDA_FOUND})
   message("CUDA was found")
   string(APPEND CMAKE_CXX_FLAGS " -DUSE_CUDA")
endif()
set(FindIPP_DIR ${CMAKE_SOURCE_DIR}/CMakeModules)
find_package(FindIPP REQUIRED)

if (MFLIB_WRAP_PYTHON)
   find_package(PythonInterp 3)
   find_package(PythonLibs 3 REQUIRED)
   message(STATUS "PYTHON_LIBRARIES = ${PYTHON_LIBRARIES}")
   message(STATUS "PYTHON_EXECUTABLE = ${PYTHON_EXECUTABLE}")
   message(STATUS "PYTHON_INCLUDE_DIRS = ${PYTHON_INCLUDE_DIRS}")
   find_package(pybind11 REQUIRED)
endif()

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${MKL_INCLUDE_DIR}
    ${IPP_INCLUDE_DIR}
    ${MPI_CXX_INCLUDE_PATH}
    ${CUDA_INCLUDE_DIRS}
)

set(SRC
    src/matchedFilterOptions.cpp
    src/matchedFilter.cpp
    src/waveformTemplate.cpp)
set(MPI_SRC
    src/mpi.cpp)
if (${CUDA_FOUND})
   set(CUDA_SRC src/cudaMatchedFilter.cpp)
endif()
#cmake -DBUILD_SHARED_LIBS=YES /path/to/source 
set(BUILD_SHARED_LIBS YES)
add_library(matchedFilter ${SRC} ${MPI_SRC} ${CUDA_SRC})

################################################################################
#                                  Unit Tests                                  #
################################################################################

add_executable(unitTests
               testing/main.cpp
               testing/waveformTemplate.cpp
               testing/matchedFilterOptions.cpp
               testing/matchedFilter.cpp)
target_link_libraries(matchedFilter ${IPP_LIBRARY} ${MPI_CXX_LIBRARIES} ${CUDA_CUFFT_LIBRARIES})
target_link_libraries(unitTests
                      PRIVATE matchedFilter ${GTEST_BOTH_LIBRARIES} ${IPP_LIBRARY} ${MPI_CXX_LIBRARIES} ${CUDA_CUFFT_LIBRARIES})
target_include_directories(unitTests PRIVATE ${GTEST_INCLUDE_DIRS})

add_test(NAME unitTests
         COMMAND unitTests)
