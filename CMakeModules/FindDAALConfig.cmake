# Already in cache, be silent
if (DAAL_INCLUDE_DIR AND DAAL_CORE_LIBRARY AND DAAL_SEQUENTIAL_LIBRARY)
   set (DAAL_FIND_QUIETLY TRUE)
endif()

#if (NOT BUILD_SHARED_LIBS)
#   set(IPPS "libipps.a")
#   set(VM "libippvm.a")
#   set(CORE "libippcore.a")
#else()
   set(CORE "daal_core")
   set(SEQUENTIAL "daal_sequential")
#endif()

find_path(DAAL_INCLUDE_DIR
          NAMES daal.h
          HINTS /opt/intel/daal/include
                $ENV{DAAL_ROOT}/include
                $ENV{DAAL_INCLUDE_DIR})
find_library(DAAL_CORE_LIBRARY
             NAMES ${CORE}
             PATHS /opt/intel/daal/lib/intel64
                   $ENV{DAAL_ROOT}/lib/intel64
                   $ENV{DAAL_ROOT}/lib/
                   $ENV{DAAL_LIB_DIR})
find_library(DAAL_SEQUENTIAL_LIBRARY
             NAMES ${SEQUENTIAL}
             PATHS /opt/intel/daal/lib/intel64
                   $ENV{DAAL_ROOT}/lib/intel64
                   $ENV{DAAL_ROOT}/lib/
                   $ENV{DAAL_LIB_DIR})

set(DAAL_LIBRARY ${DAAL_SEQUENTIAL_LIBRARY} ${DAAL_CORE_LIBRARY})

# Handle the QUIETLY and REQUIRED arguments and set MKL_FOUND to TRUE if
# all listed variables are TRUE.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(DAAL DEFAULT_MSG DAAL_LIBRARY DAAL_INCLUDE_DIR DAAL_CORE_LIBRARY DAAL_SEQUENTIAL_LIBRARY)
mark_as_advanced(DAAL_INCLUDE_DIR DAAL_LIBRARY DAAL_SEQUENTIAL_LIBRARY DAAL_CORE_LIBRARY)
