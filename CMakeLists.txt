cmake_minimum_required(VERSION 3.10)
project(mflib VERSION 1.0.0 LANGUAGES CXX)
enable_testing()

include(CheckCXXCompilerFlag)
#set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)
include(FindOpenMP QUIET)
if (OPENMP_FOUND)
   message("OpenMP found")
   string(APPEND CMAKE_CXX_FLAGS " ${OpenMP_CXX_FLAGS}")
   link_libraries(OpenMP::OpenMP_CXX)
endif()

check_cxx_source_compiles("
#include <memory>
int main () {
  std::unique_ptr<int> foo = std::make_unique<int>();
}
"
  HAVE_STD__MAKE_UNIQUE
)
check_cxx_source_compiles("
#if __has_include(<pstl/algorithm>)
#include <pstl/algorithm>
#endif
int main() {
}
"
  HAVE_SOURCE_FILE_INCLUSION
)
#if (${HAVE_STD__MAKE_UNIQUE})
#   message("Using library provided unique_ptr implementation")
#   add_definitions(matchedFilter PUBLIC HAVE_STD__MAKE_UNIQUE=1)
#   #add_definitions(-DHAVE_STD__MAKE_UNIQUE)
#endif()

find_package(GTest REQUIRED)
if (MFLIB_USE_MPI)
    find_package(MPI REQUIRED)
endif()
#find_package(CUDA
#             COMPONENTS CUDA_INCLUDE_DIRS CUDA_CUFFT_LIBRARIES)
#if (${CUDA_FOUND})
#   message("CUDA was found")
#   enable_language(CUDA)
#   string(APPEND CMAKE_CXX_FLAGS " -DUSE_CUDA")
#endif()
#set(FindIPP_DIR ${CMAKE_SOURCE_DIR}/CMakeModules)
set(FindMKL_DIR ${CMAKE_SOURCE_DIR}/CMakeModules)
set(FindDAAL_DIR ${CMAKE_SOURCE_DIR}/CMakeModules)
set(FindTBB_DIR ${CMAKE_SOURCE_DIR}/CMakeModules)
#find_package(FindIPP REQUIRED)
find_package(FindMKL REQUIRED)
find_package(FindDAAL REQUIRED)
find_package(FindTBB)

if (MFLIB_WRAP_PYTHON)
   find_package(PythonInterp 3)
   find_package(PythonLibs 3 REQUIRED)
   message(STATUS "PYTHON_LIBRARIES = ${PYTHON_LIBRARIES}")
   message(STATUS "PYTHON_EXECUTABLE = ${PYTHON_EXECUTABLE}")
   message(STATUS "PYTHON_INCLUDE_DIRS = ${PYTHON_INCLUDE_DIRS}")
   #string(APPEND CMAKE_CXX_FLAGS " -DPYBIND11_CPP17 ")
   find_package(pybind11 REQUIRED)
endif()
if (TBB_FOUND)
   message("TBB was found")
   add_definitions(-DHAVE_TBB)
endif()

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${MKL_INCLUDE_DIR}
    ${IPP_INCLUDE_DIR}
    ${MPI_CXX_INCLUDE_PATH}
    ${DAAL_INCLUDE_DIR}
    #${CUDA_INCLUDE_DIRS}
)

configure_file(${CMAKE_SOURCE_DIR}/include/mflib/version.hpp.in
               ${CMAKE_SOURCE_DIR}/include/mflib/version.hpp)
set(PUBLIC_HEADER_DIRECTORIES
    include/mflib/
    include/mflib/singleChannel/)
set(PRIVATE_HEADER_DIRECTORIES
    include/private)
set(SRC
    src/version.cpp
    src/peakFinder.cpp
    src/dbscan.cpp
    src/matchedFilterParameters.cpp
    src/matchedFilter.cpp
    src/networkStationPhase.cpp
    src/singleChannel/associator.cpp
    src/singleChannel/associatorParameters.cpp
    src/singleChannel/detectionTimeInterpolationParameters.cpp
    src/singleChannel/detectorParameters.cpp
    src/singleChannel/detector.cpp
    src/singleChannel/detection.cpp
    src/singleChannel/event.cpp
    src/singleChannel/matchedFilterParameters.cpp
    src/singleChannel/matchedFilter.cpp
    src/singleChannel/relativeMagnitude.cpp
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
target_link_libraries(matchedFilter PRIVATE ${MKL_LIBRARY} ${DAAL_LIBRARY} ${IPP_LIBRARY}
                      ${MPI_CXX_LIBRARIES} ${TBB_LIBRARY})
#                      ${CUDA_CUFFT_LIBRARIES})# OpenMP::OpenMP_CXX)
set_target_properties(matchedFilter PROPERTIES
                      CXX_STANDARD 17
                      CXX_STANDARD_REQUIRED YES
                      CXX_EXTENSIONS NO)
target_include_directories(matchedFilter PUBLIC ${PUBLIC_HEADER_FILES})
target_include_directories(matchedFilter PRIVATE ${PRIVATE_HEADER_FILES})
if(OpenMP_CXX_FOUND)
   target_link_libraries(matchedFilter PUBLIC OpenMP::OpenMP_CXX)
   #set(OPENMP_CXX_LIB OpenMP::OpenMP_CXX)
   target_include_directories(matchedFilter PRIVATE OpenMP::OpenMP_CXX)
endif()
if (${HAVE_STD__MAKE_UNIQUE})
   message("Using std::unique_ptr implementation")
   target_compile_definitions(matchedFilter PUBLIC HAVE_STD__MAKE_UNIQUE=1)
   #add_definitions(-DHAVE_STD__MAKE_UNIQUE)
endif()
if (${HAVE_SOURCE_FILE_INCLUSION})
   message("Compiler supports source file inclusion")
   target_compile_definitions(matchedFilter PUBLIC HAVE_SOURCE_FILE_INCLUSION=1)
endif()
if (MFLIB_WRAP_PYTHON)
   add_library(pymflib MODULE
               python/waveformTemplate.cpp
               python/networkStationPhase.cpp
               python/singleChannel/associator.cpp
               python/singleChannel/associatorParameters.cpp
               python/singleChannel/detection.cpp 
               python/singleChannel/detector.cpp
               python/singleChannel/detectorParameters.cpp
               python/singleChannel/matchedFilterParameters.cpp
               python/singleChannel/matchedFilter.cpp
               python/singleChannel/relativeMagnitude.cpp
               python/multiChannelMatchedFilterParameters.cpp
               python/multiChannelMatchedFilter.cpp
               python/pymflib.cpp)
   target_link_libraries(pymflib PRIVATE
                         pybind11::module matchedFilter ${MKL_LIBRARY} ${MKL_EXTRA} ${DAAL_LIBRARY} ${TBB_LIBRARY})# OpenMP::OpenMP_CXX)
   set_property(TARGET pymflib PROPERTY CXX_STD_14)
   set_property(TARGET pymflib PROPERTY PREFIX "")
   message("hey" ${CMAKE_SOURCE_DIR}/python/include)
   target_include_directories(pymflib PRIVATE pybind11 ${CMAKE_SOURCE_DIR}/python/include)
endif()

################################################################################
#                                  Unit Tests                                  #
################################################################################

add_executable(unitTests
               testing/main.cpp
               testing/networkStationPhase.cpp
               testing/waveformTemplate.cpp
               testing/matchedFilterParameters.cpp
               testing/matchedFilter.cpp
               testing/singleChannel.cpp
               testing/relativeMagnitude.cpp
               testing/peakFinder.cpp
               testing/dbscan.cpp
               testing/detection.cpp
               testing/detector.cpp
               testing/associator.cpp
               testing/utilities.cpp)
target_link_libraries(unitTests
                      PRIVATE matchedFilter ${GTEST_BOTH_LIBRARIES} ${MKL_LIBRARY} ${IPP_LIBRARY}
                      ${MPI_CXX_LIBRARIES} ${CUDA_CUFFT_LIBRARIES} ${DAAL_LIBRARY} ${TBB_LIBRARY})# OpenMP::OpenMP_CXX)
target_include_directories(unitTests PRIVATE ${GTEST_INCLUDE_DIRS})
set_target_properties(unitTests  PROPERTIES
                      CXX_STANDARD 17
                      CXX_STANDARD_REQUIRED YES 
                      CXX_EXTENSIONS NO)

file(COPY ${CMAKE_SOURCE_DIR}/python/unit_tests.py DESTINATION .)
add_test(NAME unitTests
         COMMAND unitTests)

###############################################################################
#                                 Installation                                #
###############################################################################
include(GNUInstallDirs)
if (MFLIB_WRAP_PYTHON)
   install(TARGETS matchedFilter pymflib
           RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
           LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
           ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
           PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
else()
   install(TARGETS matchedFilter
           RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
           LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
           ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
           PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
endif()

##########################################################################################
#                                     CPACK Packaging                                    #
##########################################################################################
set(CPACK_PACKAGE_NAME "mflib")
set(CPACK_PACKAGE_VENDOR "The Internet")
set(CPACK_PACKAGE_CONTACT "bbaker@seis.utah.edu")
set(CPACK_PACKAGE_LICENSE "MIT")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "A toolkit for template matching seismic signals")
set(CPACK_PACKAGE_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${PROJECT_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${PROJECT_VERSION_PATCH})
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
