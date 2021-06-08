cmake_minimum_required(VERSION 3.13)

# Project
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
project(libmolpp VERSION 0.1.0)

# Options
option(USE_SYSTEM_EIGEN "Use system's Eigen installation." OFF)
option(BUILD_SHARED_LIBS "Build shared libraries" OFF)
option(BUILD_TESTING "Build library tests." OFF)

# Third-party
add_subdirectory("thirdparty/molfile_plugin")
include(FetchContent)
if(USE_SYSTEM_EIGEN)
    find_package (Eigen3 3.4 REQUIRED NO_MODULE)
else()
    FetchContent_Declare(
        Eigen3
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen
        GIT_TAG 3.4-rc1
    )
    FetchContent_MakeAvailable(Eigen3)
endif()

# The library
add_library(molpp)

# Source files
add_subdirectory("src")

# Include dirs
target_include_directories(molpp
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
    PRIVATE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
)

# Linkings
target_link_libraries(molpp
    PUBLIC
        Eigen3::Eigen
    PRIVATE
        molfile_plugin
)

# Compiler options
target_compile_features(molpp
    PUBLIC
        cxx_std_17
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
if(BUILD_TESTING)
    add_subdirectory(tests)
endif()
