# Already in cache, be silent
if (TBB_LIBRARY)
   set (TBB_FIND_QUIETLY TRUE)
endif()

#if (NOT BUILD_SHARED_LIBS)
#   set(TBB "libtbb.a")
#else()
   set(TBB "tbb")
#endif()

find_library(TBB_LIBRARY
             NAMES ${TBB}
             PATHS /opt/intel/tbb/lib
                   /opt/intel/tbb/lib/intel64/gcc4.7
                   $ENV{TBB_ROOT}/lib/
                   $ENV{TBB_ROOT}/lib/intel64/gcc4.7
                   $ENV{TBB_LIB_DIR})

# Handle the QUIETLY and REQUIRED arguments and set MKL_FOUND to TRUE if
# all listed variables are TRUE.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(TBB DEFAULT_MSG TBB_LIBRARY)
mark_as_advanced(TBB_LIBRARY)
