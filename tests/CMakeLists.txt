find_package(GTest 1.12.1 REQUIRED)

enable_testing()

add_executable(tests
    test_Atom.cpp
    test_Property.cpp
    test_PropertyTrajectory.cpp
    test_PropertyContainer.cpp
    test_MolData.cpp
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
    dssp.cpp
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

target_link_libraries(tests
    GTest::gtest
    GTest::gtest_main
    GTest::gmock_main
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
