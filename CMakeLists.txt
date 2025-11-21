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
project(
    oncvpsp
    VERSION 4.0.1
    DESCRIPTION "Optimized norm-conserving Vanderbilt pseudopotentials"
    LANGUAGES Fortran C)

##########################################################
# Options
##########################################################
option(ONCVPSP_WITH_Libxc "Enable Libxc support" ON)
option(ONCVPSP_WITH_TOML "Enable TOML input support" ON)
option(ONCVPSP_WITH_HDF5 "Enable HDF5 output support" OFF)

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
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib
    CACHE
    PATH "Single output directory for building all libraries."
)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib
    CACHE
    PATH "Single output directory for building all libraries."
)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin
    CACHE
    PATH "Single output directory for building all executables."
)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules
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
if(DEFINED LAPACK_ROOT)
    find_package(LAPACK REQUIRED HINTS ${LAPACK_ROOT})
else()
    find_package(LAPACK REQUIRED)
endif()
message(STATUS "Found LAPACK version ${LAPACK_VERSION}")
message(STATUS "  LAPACK libraries: ${LAPACK_LIBRARIES}")
list(APPEND ONCVPSP_DEPENDENCY_LIBRARIES ${LAPACK_LIBRARIES})

###########################################################
# LIBXC
# Prefer using pkg-config; otherwise, try using CMake's FindLibxc module
###########################################################
# Try to find Libxc via pkg-config
set(Libxc_FOUND FALSE)
if(ONCVPSP_WITH_Libxc)
    if((PkgConfig_FOUND) AND (NOT DEFINED Libxc_ROOT))
        message(STATUS "Looking for Libxc via pkg-config")
        # C library
        pkg_check_modules(PkgConfig_Libxc QUIET IMPORTED_TARGET libxc)
        if(PkgConfig_Libxc_FOUND)
            message(STATUS "Found Libxc version ${PkgConfig_Libxc_VERSION} via pkg-config")
        endif()
        # Fortran library
        if(PkgConfig_Libxc_VERSION VERSION_LESS 5)
            # libxcf90 for Libxc major version 3 and 4
            pkg_check_modules(PkgConfig_LibxcFortran QUIET IMPORTED_TARGET libxcf90)
        else()
            # libxcf03 for Libxc major version 5 and above
            pkg_check_modules(PkgConfig_LibxcFortran QUIET IMPORTED_TARGET libxcf03)
        endif()
        # Fortran interface check
        if(NOT PkgConfig_LibxcFortran_FOUND)
            message(STATUS "Libxc found via pkg-config lacks Fortran interface.")
        endif()
        # Success check
        if(PkgConfig_Libxc_FOUND AND PkgConfig_LibxcFortran_FOUND)
            message(STATUS "Found Libxc via pkg-config")
            set(Libxc_FOUND TRUE)
            set(Libxc_VERSION ${PkgConfig_Libxc_VERSION})
            set(Libxc_LIBRARIES PkgConfig::PkgConfig_Libxc PkgConfig::PkgConfig_LibxcFortran)
            set(Libxc_INCLUDE_DIRS ${PkgConfig_Libxc_INCLUDE_DIRS} ${PkgConfig_LibxcFortran_INCLUDE_DIRS})
            message(STATUS "  Libxc libraries: ${PkgConfig_Libxc_LIBRARIES}")
        else()
            message(STATUS "Libxc not found via pkg-config")
        endif()
    endif()

    # If not found via pkg-config, try CMake's FindLibxc module
    if(NOT Libxc_FOUND)
        message(STATUS "Looking for Libxc via CMake's find_package")
        find_package(Libxc)
        if(Libxc_FOUND)
            message(STATUS "Found Libxc version ${Libxc_VERSION} via CMake's find_package")
            # Fortran interface check
            if(Libxc_VERSION VERSION_LESS 5)
                if(NOT TARGET Libxc::xcf90)
                    message(STATUS "Libxc found via CMake lacks Fortran90 interface.")
                endif()
            elseif(NOT TARGET Libxc::xcf03)
                message(FATAL_ERROR "Libxc found via CMake lacks Fortran03 interface.")
            endif()
        endif()
        if(Libxc_FOUND)
            if(Libxc_VERSION VERSION_LESS 5)
                set(Libxc_LIBRARIES Libxc::xcf90 Libxc::xc)
            else()
                set(Libxc_LIBRARIES Libxc::xcf03 Libxc::xc)
            endif()
            set(Libxc_INCLUDE_DIRS ${Libxc_INCLUDE_DIRS})
        endif()
    endif()

    if(NOT Libxc_FOUND)
        message(FATAL_ERROR "Libxc not found.")
    elseif(Libxc_VERSION VERSION_LESS 3)
        message(FATAL_ERROR "Libxc version ${Libxc_VERSION} is too old (<3.0.0).")
    elseif(Libxc_VERSION VERSION_EQUAL 5.0.0)
        message(FATAL_ERROR "Libxc version 5.0.0 is affected by serious bugs.")
    else()
        message(STATUS "Libxc version ${Libxc_VERSION}:")
        list(APPEND ONCVPSP_DEPENDENCY_LIBRARIES ${Libxc_LIBRARIES})
        message(STATUS "  Libxc libraries: ${Libxc_LIBRARIES}")
        list(APPEND ONCVPSP_DEPENDENCY_INCLUDE_DIRS ${Libxc_INCLUDE_DIRS})
        message(STATUS "  Libxc include directories: ${Libxc_INCLUDE_DIRS}")
    endif()
endif()

###########################################################
# TOML-F (via CPM)
###########################################################
if(ONCVPSP_WITH_TOML)
    CPMAddPackage(
        URI "gh:toml-f/toml-f@0.4.2"
        OPTIONS "BUILD_TESTING=OFF"
    )
    list(APPEND ONCVPSP_DEPENDENCY_LIBRARIES toml-f-lib)
endif()

###########################################################
# HDF5
###########################################################
set(HDF5_FOUND FALSE)
if (ONCVPSP_WITH_HDF5)
    message(STATUS "Looking for HDF5")
    find_package(HDF5 COMPONENTS Fortran HL)
    if(HDF5_USE_STATIC_LIBRARIES)
        find_package(ZLIB)
    endif()
    if (HDF5_FOUND)
        message(STATUS "Found HDF5 version ${HDF5_VERSION}")
        message(STATUS "  HDF5 include dirs: ${HDF5_INCLUDE_DIRS}")
        message(STATUS "  HDF5 libraries: ${HDF5_LIBRARIES}")
        message(STATUS "  HDF5 high-level libraries: ${HDF5_HL_LIBRARIES}")
        set(HAVE_HDF5 TRUE CACHE BOOL "Whether HDF5 is available")
        list(APPEND ONCVPSP_DEPENDENCY_LIBRARIES ${HDF5_HL_LIBRARIES})
        list(APPEND COMMON_INCLUONCVPSP_DEPENDENCY_INCLUDE_DIRSDE_DIRECTORIES ${HDF5_INCLUDE_DIRS})
    endif()
endif()

###########################################################
# Subdirectories
###########################################################
add_subdirectory(src)

###########################################################
# Installation
###########################################################
include(GNUInstallDirs)
install(TARGETS oncvpsp oncvpspnr oncvpspsr oncvpspr
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
