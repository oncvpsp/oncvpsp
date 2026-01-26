cmake_minimum_required(VERSION 3.15)

# Include Guard
set(CONFIGURED_ONCE TRUE CACHE INTERNAL
    "A flag showing that CMake has configured at least once.")

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
    if(NOT EXISTS ${LAPACK_ROOT})
        message(WARNING "LAPACK_ROOT ${LAPACK_ROOT} does not exist.")
    endif()
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
        # Fortran library
        if(PkgConfig_Libxc_VERSION VERSION_LESS 5)
            # libxcf90 for Libxc major version 3 and 4
            pkg_check_modules(PkgConfig_LibxcFortran QUIET IMPORTED_TARGET libxcf90)
        else()
            # libxcf03 for Libxc major version 5 and above
            pkg_check_modules(PkgConfig_LibxcFortran QUIET IMPORTED_TARGET libxcf03)
        endif()
        # Success check
        if(PkgConfig_Libxc_FOUND AND PkgConfig_LibxcFortran_FOUND)
            set(Libxc_FOUND TRUE)
            set(Libxc_VERSION ${PkgConfig_Libxc_VERSION})
            set(Libxc_LIBRARIES PkgConfig::PkgConfig_Libxc PkgConfig::PkgConfig_LibxcFortran)
            set(Libxc_INCLUDE_DIRS ${PkgConfig_Libxc_INCLUDE_DIRS} ${PkgConfig_LibxcFortran_INCLUDE_DIRS})
            message(STATUS "Looking for Libxc via pkg-config - found")
        else()
            message(STATUS "Looking for Libxc via pkg-config - not found")
        endif()
    endif()
    # If not found via pkg-config, try CMake's FindLibxc module
    if(NOT Libxc_FOUND)
        message(STATUS "Looking for Libxc via CMake's find_package")
        if ((DEFINED Libxc_ROOT) AND NOT (EXISTS ${Libxc_ROOT}))
            message(WARNING "Libxc_ROOT ${Libxc_ROOT} does not exist.")
        endif()
        find_package(Libxc)
        if(Libxc_FOUND)
            # Fortran interface check
            if(Libxc_VERSION VERSION_LESS 5)
                if(NOT TARGET Libxc::xcf90)
                    message(FATAL_ERROR "Libxc version ${Libxc_VERSION} found via CMake lacks the required Fortran90 interface.")
                endif()
            elseif(NOT TARGET Libxc::xcf03)
                message(FATAL_ERROR "Libxc version ${Libxc_VERSION} found via CMake lacks the required Fortran03 interface.")
            endif()
        endif()
        if(Libxc_FOUND)
            message(STATUS "Looking for Libxc via CMake's find_package - found")
            if(Libxc_VERSION VERSION_LESS 5)
                set(Libxc_LIBRARIES Libxc::xcf90 Libxc::xc)
            else()
                set(Libxc_LIBRARIES Libxc::xcf03 Libxc::xc)
            endif()
            set(Libxc_INCLUDE_DIRS ${Libxc_INCLUDE_DIRS})
        else()
            message(STATUS "Looking for Libxc via CMake's find_package - not found")
        endif()
    endif()

    if(NOT Libxc_FOUND)
        message(FATAL_ERROR "Libxc not found.")
    elseif(Libxc_VERSION VERSION_LESS 3)
        message(FATAL_ERROR "Libxc version ${Libxc_VERSION} is too old (<3.0.0).")
    elseif(Libxc_VERSION VERSION_EQUAL 5.0.0)
        message(FATAL_ERROR "Libxc version 5.0.0 is affected by serious bugs.")
    else()
        message(STATUS "Found Libxc version ${Libxc_VERSION}")
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
    set(TEMP_BUILD_TESTING ${BUILD_TESTING})
    set(BUILD_TESTING OFF)
    CPMAddPackage(
        URI "gh:toml-f/toml-f@0.4.3"
        OPTIONS "TOML_F_BUILD_TESTING=OFF CMAKE_WARN_DEPRECATED=OFF"
    )
    list(APPEND ONCVPSP_DEPENDENCY_LIBRARIES toml-f-lib)
    set(BUILD_TESTING ${TEMP_BUILD_TESTING})
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
