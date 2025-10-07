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
# PkgConfig (used preferentially for Libxc)
###########################################################
find_package(PkgConfig)

###########################################################
# LAPACK
###########################################################
find_package(LAPACK)
if (NOT LAPACK_FOUND)
    message(STATUS "LAPACK not found. Trying with hints.")
    find_package(LAPACK REQUIRED HINTS ${LAPACK_ROOT})
endif()
list(APPEND COMMON_LIBRARIES ${LAPACK_LIBRARIES})

###########################################################
# LIBXC
# Prefer using pkg-config
# Otherwise, try using CMake's FindLibxc module
###########################################################
include(CMakePrintHelpers)
if ((PkgConfig_FOUND) AND (NOT DEFINED Libxc_ROOT))
    pkg_check_modules(ONCVPSP_LIBXCF90 QUIET IMPORTED_TARGET libxcf90)
    pkg_check_modules(ONCVPSP_LIBXCF03 QUIET IMPORTED_TARGET libxcf03)
    if (ONCVPSP_LIBXCF03_FOUND AND ONCVPSP_LIBXCF90_FOUND)
        message(STATUS "Found Libxc via pkg-config")
        set(Libxc_FOUND TRUE)
        set(Libxc_VERSION ${ONCVPSP_LIBXCF03_VERSION})
        set(LIBXC_LIBRARIES PkgConfig::ONCVPSP_LIBXCF90 PkgConfig::ONCVPSP_LIBXCF03)
        set(LIBXC_INCLUDE_DIRS ${ONCVPSP_LIBXCF03_INCLUDE_DIRS})
    else()
        message(STATUS "Libxc not found via pkg-config")
        set(Libxc_FOUND FALSE)
    endif()
endif()

if (NOT Libxc_FOUND)
    find_package(Libxc COMPONENTS Fortran)
    message(STATUS "Found Libxc via CMake's find_package")
    if (Libxc_FOUND)
        set(LIBXC_LIBRARIES Libxc::xcf90 Libxc::xcf03)
        set(LIBXC_INCLUDE_DIRS ${Libxc_INCLUDE_DIRS})
    endif()
endif()

if (Libxc_FOUND)
    message(STATUS "Libxc version ${Libxc_VERSION}:")
    if (${Libxc_VERSION} VERSION_GREATER_EQUAL 5 AND ${Libxc_VERSION} VERSION_LESS 7)
        list(APPEND COMMON_LIBRARIES ${LIBXC_LIBRARIES})
        message(STATUS "  Libxc libraries: ${LIBXC_LIBRARIES}")
        list(APPEND COMMON_INCLUDE_DIRECTORIES ${LIBXC_INCLUDE_DIRS})
        message(STATUS "  Libxc include directories: ${LIBXC_INCLUDE_DIRS}")
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
