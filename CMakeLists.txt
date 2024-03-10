cmake_minimum_required(VERSION 3.13)

# Options
option(BUILD_SHARED_LIBS "Build shared libraries" OFF)
option(MOLPP_BUILD_TESTING "Build library tests." OFF)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(THREADS_PREFER_PTHREAD_FLAG ON)

project(libmolpp
    VERSION
        0.1.0
    LANGUAGES
        CXX C
)

# Dependencies
find_package(Threads REQUIRED)
find_package(Eigen3 3.4 REQUIRED NO_MODULE)
find_package(peglib REQUIRED)

# The library
add_library(molpp)
add_library(Molpp::molpp ALIAS molpp)

# Source files
add_subdirectory("include")
add_subdirectory("src")

# Third-party
add_subdirectory("thirdparty")

# Molfile_plugin definitions. Must be set in the same
# file as the molpp declaration.
foreach(FILE ${MOLFILE_SOURCES})
    get_filename_component(NAME ${FILE} NAME_WLE)
    set_property(SOURCE ${FILE} PROPERTY
        COMPILE_DEFINITIONS
            VMDPLUGIN=${NAME}
    )
endforeach(FILE)

target_link_libraries(molpp
    PUBLIC
        Eigen3::Eigen
        peglib
    PRIVATE
        Threads::Threads
)

if (MSVC)
    add_compile_options(/W4 /WX)
else()
    add_compile_options(-Wall -Wextra -pedantic -Werror)
endif()

# Testing
if(MOLPP_BUILD_TESTING)
    find_package(GTest 1.12.1 REQUIRED)
    add_subdirectory(tests)
endif()

# Installation
install(TARGETS molpp
    EXPORT
        MolppTargets
    FILE_SET molpp_headers DESTINATION
        include
    LIBRARY
        DESTINATION lib
    ARCHIVE
        DESTINATION lib
)

# Exporting targets
include(CMakePackageConfigHelpers)

SET(EXPORT_DIR "lib/cmake/Molpp")

export(
    TARGETS
        molpp
    NAMESPACE
        Molpp::
    FILE
        MolppTargets.cmake
)

configure_package_config_file(cmake/MolppConfig.cmake.in
    MolppConfig.cmake
    INSTALL_DESTINATION
        ${EXPORT_DIR}
)

write_basic_package_version_file(
    MolppConfigVersion.cmake
    VERSION
        ${PACKAGE_VERSION}
    COMPATIBILITY
        AnyNewerVersion
)

install(
    EXPORT
        MolppTargets
    FILE
        MolppTargets.cmake
    NAMESPACE
        Molpp::
    DESTINATION
        ${EXPORT_DIR}
)

install(
    FILES
        "${CMAKE_CURRENT_BINARY_DIR}/MolppConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/MolppConfigVersion.cmake"
    DESTINATION
        ${EXPORT_DIR}
)
