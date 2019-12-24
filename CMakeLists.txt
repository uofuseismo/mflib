cmake_minimum_required(VERSION 3.10)
project(mflib VERSION 0.1.3 LANGUAGES CXX)
enable_testing()

include(CheckCXXCompilerFlag)
#set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)
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
   enable_language(CUDA)
   string(APPEND CMAKE_CXX_FLAGS " -DUSE_CUDA")
endif()
set(FindIPP_DIR ${CMAKE_SOURCE_DIR}/CMakeModules)
set(FindMKL_DIR ${CMAKE_SOURCE_DIR}/CMakeModules)
find_package(FindIPP REQUIRED)
find_package(FindMKL REQUIRED)

if (MFLIB_WRAP_PYTHON)
   find_package(PythonInterp 3)
   find_package(PythonLibs 3 REQUIRED)
   message(STATUS "PYTHON_LIBRARIES = ${PYTHON_LIBRARIES}")
   message(STATUS "PYTHON_EXECUTABLE = ${PYTHON_EXECUTABLE}")
   message(STATUS "PYTHON_INCLUDE_DIRS = ${PYTHON_INCLUDE_DIRS}")
   #string(APPEND CMAKE_CXX_FLAGS " -DPYBIND11_CPP17 ")
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
#    src/shiftAndStack.cpp)
# Use MPI?
if (${MPI_FOUND})
    message("Building with MPI")
    set(SRC ${SRC}
        src/mpi.cpp)
endif()
if (${CUDA_FOUND})
   string(APPEND CUDA_NVCC_FLAGS " -ccbin ${CMAKE_CXX_COMPILER}")
   set(CUDA_SRC src/cudaMatchedFilter.cu)
endif()
#cmake -DBUILD_SHARED_LIBS=YES /path/to/source 
set(BUILD_SHARED_LIBS YES)
add_library(matchedFilter ${SRC} ${MPI_SRC} ${CUDA_SRC})
SET_PROPERTY(TARGET matchedFilter PROPERTY CXX_STANDARD 17)
#if(OpenMP_CXX_FOUND)
#   message("yes")
#   #target_link_libraries(matchedFilter PUBLIC OpenMP::OpenMP_CXX)
#   set(OPENMP_CXX_LIB OpenMP::OpenMP_CXX)
#   #message(${OPENMP_CXX_LIB})
#   target_include_directories(matchedFilter PRIVATE OpenMP::OpenMP_CXX)
#endif()
if (MFLIB_WRAP_PYTHON)
   add_library(pymflib MODULE
               python/waveformTemplate.cpp
               python/pymflib.cpp)
   target_link_libraries(pymflib PRIVATE pybind11::module matchedFilter)
   set_property(TARGET pymflib PROPERTY CXX_STD_14)
   set_property(TARGET pymflib PROPERTY PREFIX "")
   #target_include_directories(pymflib PRIVATE pybind11)
endif()

################################################################################
#                                  Unit Tests                                  #
################################################################################

add_executable(unitTests
               testing/main.cpp
               testing/waveformTemplate.cpp
               testing/matchedFilterOptions.cpp
               testing/matchedFilter.cpp)
target_link_libraries(matchedFilter ${MKL_LIBRARY} ${IPP_LIBRARY} ${MPI_CXX_LIBRARIES}
                      ${CUDA_CUFFT_LIBRARIES} OpenMP::OpenMP_CXX)
target_link_libraries(unitTests
                      PRIVATE matchedFilter ${GTEST_BOTH_LIBRARIES} ${IPP_LIBRARY} ${MPI_CXX_LIBRARIES} ${CUDA_CUFFT_LIBRARIES})
target_include_directories(unitTests PRIVATE ${GTEST_INCLUDE_DIRS})

file(COPY ${CMAKE_SOURCE_DIR}/python/unit_tests.py DESTINATION .)
add_test(NAME unitTests
         COMMAND unitTests)
