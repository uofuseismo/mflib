# Already in cache, be silent
if (IPP_INCLUDE_DIR AND IPP_IPPS_LIBRARY AND IPP_VM_LIBRARY AND IPP_CORE_LIBRARY)
   set (IPP_FIND_QUIETLY TRUE)
endif()

#if (NOT BUILD_SHARED_LIBS)
#   set(IPPS "libipps.a")
#   set(VM "libippvm.a")
#   set(CORE "libippcore.a")
#else()
   set(IPPS "ipps")
   set(VM "ippvm")
   set(CORE "ippcore")
#endif()

find_path(IPP_INCLUDE_DIR
          NAMES ipps.h
          HINTS /opt/intel/ipp/include
                $ENV{IPP_ROOT}/include
                $ENV{IPP_INC_DIR})
find_library(IPP_IPPS_LIBRARY
             NAMES ${IPPS}
             PATHS /opt/intel/ipp//lib/intel64
                   $ENV{IPP_ROOT}/lib/intel64
                   $ENV{IPP_ROOT}/lib/
                   $ENV{IPP_LIB_DIR})
find_library(IPP_VM_LIBRARY
             NAMES ${VM}
             PATHS /opt/intel/ipp/lib/intel64
                   $ENV{IPP_ROOT}/lib/intel64
                   $ENV{IPP_ROOT}/lib/
                   $ENV{IPP_LIB_DIR})
find_library(IPP_CORE_LIBRARY
             NAMES ${CORE}
             PATHS /opt/intel/ipp/lib/intel64
                   $ENV{IPP_ROOT}/lib/intel64
                   $ENV{IPP_ROOT}/lib/
                   $ENV{IPP_LIB_DIR})

set(IPP_LIBRARY ${IPP_IPPS_LIBRARY} ${IPP_VM_LIBRARY} ${IPP_CORE_LIBRARY})

# Handle the QUIETLY and REQUIRED arguments and set MKL_FOUND to TRUE if
# all listed variables are TRUE.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(IPP DEFAULT_MSG IPP_LIBRARY IPP_INCLUDE_DIR IPP_IPPS_LIBRARY IPP_VM_LIBRARY IPP_CORE_LIBRARY)
mark_as_advanced(IPP_INCLUDE_DIR IPP_LIBRARY IPP_IPPS_LIBRARY IPP_VM_LIBRARY IPP_CORE_LIBRARY)
