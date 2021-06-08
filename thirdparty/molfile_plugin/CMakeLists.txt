cmake_minimum_required(VERSION 3.13)

# Project
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
project(molfile_plugin)

# Library definitions
add_subdirectory("src")

# Include dirs
target_include_directories(molfile_plugin SYSTEM
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
    PRIVATE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
)

set_target_properties(molfile_plugin PROPERTIES
    POSITION_INDEPENDENT_CODE ${BUILD_SHARED_LIBS}
)

# Compiler options
set_target_properties(molfile_plugin PROPERTIES
    C_EXTENSIONS OFF
    CXX_EXTENSIONS OFF
)
