enable_testing()

add_executable(tests
    test_Atom.cpp
    test_MolfileReader.cpp
    test_Property.cpp
    test_PropertyTrajectory.cpp
    test_PropertyContainer.cpp
    test_MolData.cpp
    test_AtomAggregate.cpp
    test_auxiliary.cpp
    test_Graph.cpp
    test_SimpleGraph.cpp
    test_BreadthFirstTraversal.cpp
    test_ConnectedComponents.cpp
    atoms.cpp
    residues.cpp
    selections.cpp
    selectors.cpp
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
        $<TARGET_PROPERTY:molpp,INCLUDE_DIRECTORIES>
)

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

include(GoogleTest)
gtest_discover_tests(tests)
