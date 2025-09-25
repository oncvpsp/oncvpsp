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
    VERSION 4.0.1
    DESCRIPTION "Optimized norm-conserving Vanderbilt pseudopotentials"
    LANGUAGES Fortran C)

##########################################################
# Build types
##########################################################
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "Building in debug mode")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -O0 -Wall -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow")
endif (CMAKE_BUILD_TYPE STREQUAL "Debug")
if (CMAKE_BUILD_TYPE STREQUAL "Release")
    message(STATUS "Building in release mode")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -ffast-math -funroll-loops -DNDEBUG")
endif (CMAKE_BUILD_TYPE STREQUAL "Release")

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
set(COMMON_LIBRARIES "")
set(COMMON_INCLUDE_DIRECTORIES "")

###########################################################
# LAPACK
###########################################################
find_package(LAPACK REQUIRED HINTS ${LAPACK_ROOT} ${BLAS_ROOT})
if (LAPACK_FOUND)
    message(STATUS "LAPACK found.")
    message(STATUS "  LAPACK libraries: ${LAPACK_LIBRARIES}")
    message(STATUS "  LAPACK linker flags: ${LAPACK_LINKER_FLAGS}")
else ()
    message(FATAL ERROR "Failed to find LAPACK.")
endif()
list(APPEND COMMON_LIBRARIES ${LAPACK_LIBRARIES})

###########################################################
# LIBXC
###########################################################
find_package(Libxc COMPONENTS Fortran)
if (Libxc_FOUND)
    message(STATUS "Libxc version ${Libxc_VERSION} found.")
    if (${Libxc_VERSION} VERSION_GREATER_EQUAL 5 AND ${Libxc_VERSION} VERSION_LESS 7)
        set(LIBXC_LIBRARIES Libxc::xcf03)
        list(APPEND COMMON_LIBRARIES ${LIBXC_LIBRARIES})
        message(STATUS "  Libxc libraries: ${LIBXC_LIBRARIES}")
    else()
        message(STATUS "Only Libxc major versions 5 and 6 are supported. Setting Libxc_FOUND to FALSE.")
        set(Libxc_FOUND FALSE)
    endif()
endif()

###########################################################
# CMake Configuration
###########################################################
# Save compile commands to `compile_commands.json` file
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
link_directories(${CMAKE_LIBRARY_PATH} ${COMMON_INCLUDE_DIRECTORIES})
include_directories(${CMAKE_INCLUDE_PATH} ${COMMON_INCLUDE_DIRECTORIES})

###########################################################
# Subdirectories
###########################################################
add_subdirectory(src)
