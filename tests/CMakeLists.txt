enable_testing()

add_executable(tests
    test_Atom.cpp
    test_Residue.cpp
    test_MolfileReader.cpp
    test_MolData.cpp
    test_MolecularEntity.cpp
    test_Graph.cpp
    test_SimpleGraph.cpp
    test_BreadthFirstTraversal.cpp
    test_ConnectedComponents.cpp
    test_ResidueDetect.cpp
    test_ResidueData.cpp
    test_SelIndices.cpp
    test_Sel.cpp
    test_conversions.cpp
    test_ChainData.cpp
    test_AtomSel.cpp
    test_ResidueSel.cpp
    test_Topology.cpp
    atoms.cpp
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
