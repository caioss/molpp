cmake_minimum_required(VERSION 3.13)

# Project
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
project(libmolpp VERSION 0.1.0)

# Options
option(MOLPP_USE_SYSTEM_EIGEN "Use system's Eigen installation." OFF)
option(BUILD_SHARED_LIBS "Build shared libraries" OFF)
option(MOLPP_BUILD_TESTING "Build library tests." OFF)

# The library
add_library(molpp)
add_library(Molpp::Molpp ALIAS molpp)

# Source files
add_subdirectory("src")

# Third-party
add_subdirectory("thirdparty")

include(FetchContent)
if(MOLPP_USE_SYSTEM_EIGEN)
    find_package(Eigen3 3.4 REQUIRED NO_MODULE)
else()
    FetchContent_Declare(
        Eigen3
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen
        GIT_TAG 3.4.0
    )
    FetchContent_MakeAvailable(Eigen3)
endif()

# Molfile_plugin definitions. Must be set in the same
# file as the molpp declaration.
foreach(FILE ${MOLFILE_SOURCES})
    get_filename_component(NAME ${FILE} NAME_WLE)
    set_property(SOURCE ${FILE} PROPERTY
        COMPILE_DEFINITIONS
            VMDPLUGIN=${NAME}
    )
endforeach(FILE)

# Include dirs
target_include_directories(molpp
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
    PRIVATE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
)

# Linkings
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads)

target_link_libraries(molpp
    PUBLIC
        Eigen3::Eigen
    PRIVATE
        Threads::Threads
)

# Compiler options
target_compile_features(molpp
    PUBLIC
        cxx_std_20
)
set_target_properties(molpp PROPERTIES
    CXX_EXTENSIONS OFF
)
set(gcc_like_cxx "$<COMPILE_LANG_AND_ID:CXX,ARMClang,AppleClang,Clang,GNU>")
set(msvc_cxx "$<COMPILE_LANG_AND_ID:CXX,MSVC>")
target_compile_options(molpp
    PRIVATE
        "$<${gcc_like_cxx}:$<BUILD_INTERFACE:-Wall;-Wextra;-Wshadow;-pedantic>>"
        "$<${msvc_cxx}:$<BUILD_INTERFACE:-W3>>"
)

# Testing
if(MOLPP_BUILD_TESTING)
    add_subdirectory(tests)
endif()

# Exporting files
export(TARGETS molpp NAMESPACE Molpp:: FILE MolppTargets.cmake)
configure_file(cmake/MolppConfig.cmake.in MolppConfig.cmake @ONLY)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    MolppConfigVersion.cmake
    VERSION
        ${PACKAGE_VERSION}
    COMPATIBILITY
        AnyNewerVersion
)

# Installation
install(TARGETS molpp
    EXPORT
        MolppTargets
    LIBRARY
        DESTINATION lib
    ARCHIVE
        DESTINATION lib
)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/
    DESTINATION include
)

install(
    EXPORT
        MolppTargets
    FILE
        MolppTargets.cmake
    NAMESPACE
        Molpp::
    DESTINATION
        lib/cmake/Molpp
)

install(
    FILES
        "${CMAKE_CURRENT_BINARY_DIR}/MolppConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/MolppConfigVersion.cmake"
    DESTINATION
        lib/cmake/Molpp
)
