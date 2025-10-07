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

###########################################################
# Include CMake modules
###########################################################
include(cmake/CPM.cmake)

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
if (DEFINED LAPACK_ROOT)
    find_package(LAPACK REQUIRED HINTS ${LAPACK_ROOT})
else()
    find_package(LAPACK REQUIRED)
endif()
message(STATUS "Found LAPACK version ${LAPACK_VERSION}")
message(STATUS "  LAPACK libraries: ${LAPACK_LIBRARIES}")
list(APPEND COMMON_LIBRARIES ${LAPACK_LIBRARIES})

###########################################################
# LIBXC
# Prefer using pkg-config; otherwise, try using CMake's FindLibxc module
# N.B. Libxc v5.0.0 had a breaking change to the Fortran interface where
# the f03 module was renamed (and overwrote an existing) f90 module, so
# it needs to be special-cased.
###########################################################

# Try to find Libxc via pkg-config
set(Libxc_FOUND FALSE)
if ((PkgConfig_FOUND) AND (NOT DEFINED Libxc_ROOT))
    message(STATUS "Looking for Libxc via pkg-config")
    # C library
    pkg_check_modules(ONCVPSP_LIBXC QUIET IMPORTED_TARGET libxc)
    if (ONCVPSP_LIBXC_FOUND)
        message(STATUS "Found Libxc version ${ONCVPSP_LIBXC_VERSION} via pkg-config")
    endif()
    # Fortran interface
    pkg_check_modules(ONCVPSP_LIBXCF QUIET IMPORTED_TARGET libxcf03)
    if (ONCVPSP_LIBXCF_FOUND)
        message(STATUS "Found Libxc Fortran version ${ONCVPSP_LIBXCF_VERSION} via pkg-config")
    endif()
    # Version check
    if (ONCVPSP_LIBXC_VERSION VERSION_LESS 3)
        message(STATUS "Libxc version ${ONCVPSP_LIBXC_VERSION} found via pkg-config is too old (<3.0.0). Setting Libxc_FOUND to FALSE.")
        set(Libxc_FOUND FALSE)
    elseif (ONCVPSP_LIBXC_VERSION VERSION_EQUAL 5.0.0)
        message(STATUS "Libxc version 5.0.0 found via pkg-config is affected by serious bugs. Setting Libxc_FOUND to FALSE.")
        set(Libxc_FOUND FALSE)
    endif()
    # Success check
    if (ONCVPSP_LIBXC_FOUND AND ONCVPSP_LIBXCF_FOUND)
        message(STATUS "Found Libxc via pkg-config")
        set(Libxc_FOUND TRUE)
        set(Libxc_VERSION ${ONCVPSP_LIBXCF_VERSION})
        set(LIBXC_LIBRARIES PkgConfig::ONCVPSP_LIBXC PkgConfig::ONCVPSP_LIBXCF)
        set(LIBXC_INCLUDE_DIRS ${ONCVPSP_LIBXCF_INCLUDE_DIRS})
    else()
        message(STATUS "Libxc not found via pkg-config")
    endif()
endif()

# If not found via pkg-config, try CMake's FindLibxc module
if (NOT Libxc_FOUND)
    message(STATUS "Looking for Libxc via CMake's find_package")
    find_package(Libxc)
    if (Libxc_FOUND)
        message(STATUS "Found Libxc version ${Libxc_VERSION} via CMake's find_package")
        # Version check
        if (Libxc_VERSION VERSION_LESS 3)
            message(STATUS "Libxc version ${Libxc_VERSION} found via CMake is too old. Setting Libxc_FOUND to FALSE.")
            set(Libxc_FOUND FALSE)
        elseif (Libxc_VERSION VERSION_EQUAL 5.0.0)
            message(STATUS "Libxc version 5.0.0 found via pkg-config is affected by serious bugs. Setting Libxc_FOUND to FALSE.")
            set(Libxc_FOUND FALSE)
        endif()
        # Fortran interface check
        if (NOT TARGET Libxc::xcf03)
            message(STATUS "Libxc version ${Libxc_VERSION} found via CMake lacks Fortran03 interface. Setting Libxc_FOUND to FALSE.")
            set(Libxc_FOUND FALSE)
        elseif (Libxc_VERSION VERSION_LESS 7)
            if (NOT TARGET Libxc::xcf90)
                message(STATUS "Libxc version ${Libxc_VERSION} found via CMake lacks Fortran90 interface. Setting Libxc_FOUND to FALSE.")
                set(Libxc_FOUND FALSE)
            endif()
        endif()
    endif()
    if (Libxc_FOUND)
        if (Libxc_VERSION VERSION_GREATER_EQUAL 7)
            set(LIBXC_LIBRARIES Libxc::xcf03 Libxc::xc)
        else()
            set(LIBXC_LIBRARIES Libxc::xcf03 Libxc::xcf90 Libxc::xc)
        endif()
        set(LIBXC_INCLUDE_DIRS ${Libxc_INCLUDE_DIRS})
    endif()
endif()

if (Libxc_FOUND)
    message(STATUS "Libxc version ${Libxc_VERSION}:")
    list(APPEND COMMON_LIBRARIES ${LIBXC_LIBRARIES})
    message(STATUS "  Libxc libraries: ${LIBXC_LIBRARIES}")
    list(APPEND COMMON_INCLUDE_DIRECTORIES ${LIBXC_INCLUDE_DIRS})
    message(STATUS "  Libxc include directories: ${LIBXC_INCLUDE_DIRS}")
endif()

###########################################################
# TOML-F (via CPM)
###########################################################
CPMAddPackage(
    URI "gh:toml-f/toml-f@0.4.2"
    OPTIONS "BUILD_TESTING=OFF"
)
list(APPEND COMMON_LIBRARIES toml-f-lib)

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
