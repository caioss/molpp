include(FetchContent)
FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest
    GIT_TAG release-1.11.0
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()

add_executable(tests
    atoms.cpp
    residues.cpp
    atomaggregates.cpp
    selections.cpp
    selectors.cpp
    readers.cpp
    system.cpp
    bonds.cpp
    tools.cpp
    tables.cpp
    files.cpp
    auxiliary.cpp
)

target_include_directories(tests
    PRIVATE
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/thirdparty/cpp-peglib>
)

# Needed by cpp-peblib
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads)

target_link_libraries(
    tests
    gtest_main
    gmock_main
    molpp
    Threads::Threads
)

target_compile_features(tests
    PUBLIC
        cxx_std_17
)

set_target_properties(tests PROPERTIES
    CXX_EXTENSIONS OFF
)

include(GoogleTest)
gtest_discover_tests(tests)
