cmake_minimum_required(VERSION 3.15)

# Include Guard
set(CONFIGURED_ONCE TRUE CACHE INTERNAL
    "A flag showing, that CMake has configured at least once.")

# Prevent in-source builds
if(${CMAKE_BINARY_DIR} STREQUAL ${CMAKE_SOURCE_DIR})
    message(FATAL ERROR "Source folder cannot be safely used as build folder!")
endif()

##########################################################
# oncvpsp
##########################################################
project(oncvpsp
    VERSION 2.1.1
    DESCRIPTION "Optimized norm-conserving Vanderbilt pseudopotentials"
    LANGUAGES Fortran C)

###########################################################
# Include CMake modules
###########################################################
include(cmake/CPM.cmake)

##########################################################
# Build types
##########################################################
set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -fPIC -Wall -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -fPIC -ffast-math -funroll-loops -DNDEBUG")

##########################################################
# Define the paths for static libraries and executables
##########################################################
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
    CACHE
    PATH "Single output directory for building all libraries."
)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
    CACHE
    PATH "Single output directory for building all libraries."
)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
    CACHE
    PATH "Single output directory for building all executables."
)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules
    CACHE
    PATH "Single output directory for building all Fortran modules."
)

###########################################################
# Common variables
###########################################################
set(ONCVPSP_DEPENDENCY_LIBRARIES "")
set(ONCVPSP_DEPENDENCY_INCLUDE_DIRS "")

###########################################################
# PkgConfig (used preferentially for Libxc)
###########################################################
find_package(PkgConfig)

###########################################################
# LAPACK
###########################################################
if (DEFINED LAPACK_ROOT)
    find_package(LAPACK REQUIRED HINTS ${LAPACK_ROOT})
else()
    find_package(LAPACK REQUIRED)
endif()
message(STATUS "Found LAPACK version ${LAPACK_VERSION}")
message(STATUS "  LAPACK libraries: ${LAPACK_LIBRARIES}")
list(APPEND ONCVPSP_DEPENDENCY_LIBRARIES ${LAPACK_LIBRARIES})

###########################################################
# Subdirectories
###########################################################
add_subdirectory(src)

###########################################################
# Installation
###########################################################
include(GNUInstallDirs)
install(TARGETS oncvpsp oncvpspnr oncvpspsr oncvpspfr
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
