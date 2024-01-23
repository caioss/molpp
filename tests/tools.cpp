#include "auxiliary.hpp"
#include <molpp/tools/iterators.hpp>
#include "tools/algorithms.hpp"
#include <molpp/tools/Graph.hpp>
#include "tools/math.hpp"
#include "tools/SpatialSearch.hpp"
#include <molpp/MolppCore.hpp>
#include <molpp/internal/SelIndex.hpp>
#include <molpp/internal/VectorView.hpp>
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <vector>
#include <numeric>

using namespace mol::internal;
using namespace testing;

TEST(Iterators, IteratorWrapper) {
    using iter_type = typename std::vector<int>::iterator;
    std::vector<int> values(2);

    class Wrapper : public IteratorWrapper<iter_type>
    {
        using base = IteratorWrapper<iter_type>;
        using base::IteratorWrapper;
    };

    Wrapper iter1(values.begin());
    Wrapper iter2(values.begin());
    Wrapper end(values.end());

    EXPECT_EQ(iter1, iter2);
    EXPECT_EQ(iter1++, iter2);
    EXPECT_TRUE(iter1 != iter2);
    EXPECT_EQ(++iter1, end);
    EXPECT_EQ(iter1 - iter2, 2);
}

TEST(Algorithms, BreadthFirstSearch) {
    using GraphInt = Graph<int, int>;
    GraphInt graph;

    for (int i = 0; i < 6; ++i)
    {
        ASSERT_TRUE(graph.add_node(i));
    }
    graph.add_edge(0, 1, 0);
    graph.add_edge(0, 2, 0);
    graph.add_edge(2, 3, 0);
    graph.add_edge(1, 4, 0);
    graph.add_edge(4, 5, 0);
    graph.add_edge(5, 3, 0);

    BreadthFirstTraversal bfs(graph);
    EXPECT_TRUE(bfs.run(0, [](int const &node) {
        return node == 3;
    }, [](auto){
        return true;
    }));
    // Node 4 may or may not be present (STL-implementation-dependent)
    EXPECT_THAT(bfs.visited(), IsSupersetOf({0, 1, 2, 3}));
    EXPECT_THAT(bfs.visited(), IsSubsetOf({0, 1, 2, 3, 4}));
    EXPECT_THAT(bfs.parent_map(), IsSupersetOf({Pair(3, 2), Pair(2, 0), Pair(1, 0)}));
    EXPECT_THAT(bfs.parent_map(), IsSubsetOf({Pair(3, 2), Pair(2, 0), Pair(4, 1), Pair(1, 0)}));

    // Search using a mask (without node 4)
    EXPECT_TRUE(bfs.run(0, [](int const &node) {
        return node == 3;
    }, [](int const &node){
        return node != 4;
    }));
    EXPECT_THAT(bfs.visited(), UnorderedElementsAre(0, 1, 2, 3));
    EXPECT_THAT(bfs.parent_map(), UnorderedElementsAre(Pair(3, 2), Pair(2, 0), Pair(1, 0)));

    // Repeat search but starting at node 4
    EXPECT_FALSE(bfs.run(4, [](int const &node) {
        return node == 3;
    }, [](int const &node){
        return node != 4;
    }));

    // Search not possible
    EXPECT_FALSE(bfs.run(0, [](int const &node) {
        return node == 3;
    }, [](int const &node){
        return node == 0 || node == 3 || node == 5;
    }));
}

TEST(Algorithms, ConnectedComponents) {
    using GraphInt = Graph<int, int>;
    GraphInt graph;

    for (int i = 0; i < 7; ++i)
    {
        ASSERT_TRUE(graph.add_node(i));
    }

    ConnectedComponents components(graph);

    // Without edges
    EXPECT_EQ(components.run([](auto){return true;}), 7);
    EXPECT_THAT(components.components(), UnorderedElementsAre(ElementsAre(0), ElementsAre(1), ElementsAre(2), ElementsAre(3), ElementsAre(4), ElementsAre(5), ElementsAre(6)));

    graph.add_edge(0, 1, 0);
    graph.add_edge(0, 2, 0);
    graph.add_edge(2, 3, 0);
    graph.add_edge(4, 5, 0);

    // With edges
    EXPECT_EQ(components.run([](auto){return true;}), 3);
    EXPECT_THAT(components.components(), UnorderedElementsAre(UnorderedElementsAre(0, 1, 2, 3), UnorderedElementsAre(4, 5), ElementsAre(6)));

    // Filter out node 2
    EXPECT_EQ(components.run([](auto const &node){return node != 2;}), 4);
    EXPECT_THAT(components.components(), UnorderedElementsAre(UnorderedElementsAre(0, 1), ElementsAre(3), UnorderedElementsAre(4, 5), ElementsAre(6)));
}

TEST(DataStructures, SpatialSearch) {
    Eigen::Matrix3Xf points(3, 10);
    points << 0.394864, -1.92212, 0.507918, -1.10244, 2.86075, -2.68154, 2.28578, 2.88865, -1.65935, 2.21034, 0.665579, 1.90012, -0.467064, -2.63234, 2.86835, -1.38607, 0.375439, 2.77146, -1.7135, -1.00928, 0.0346084, -1.89917, -2.84799, -2.49805, 2.24333, -2.44597, 0.981083, -2.16332, -2.92808, 0.216672;

    SpatialSearch<Eigen::Matrix3Xf> search(points, 3.5);

    EXPECT_THAT(search.query(3, 3.0), UnorderedElementsAre(2, 3, 5, 8));
    EXPECT_THAT(search.query(3, 5.0), UnorderedElementsAre(0, 1, 2, 3, 5, 8, 9));

    EXPECT_THAT(search.pairs(3.0), UnorderedElementsAre(
        FieldsAre(3, 2, FloatNear(7.4041, 0.0001)),
        FieldsAre(5, 3, FloatNear(4.0494, 0.0001)),
        FieldsAre(6, 0, FloatNear(4.5555, 0.0001)),
        FieldsAre(6, 4, FloatNear(8.1384, 0.0001)),
        FieldsAre(8, 2, FloatNear(6.2570, 0.0001)),
        FieldsAre(8, 3, FloatNear(1.3393, 0.0001)),
        FieldsAre(8, 5, FloatNear(1.3845, 0.0001)),
        FieldsAre(9, 0, FloatNear(6.1342, 0.0001)),
        FieldsAre(9, 6, FloatNear(2.5074, 0.0001))));

    EXPECT_THAT(search.pairs(4.0), UnorderedElementsAre(
        FieldsAre(1, 0, FloatNear(10.6320, 0.0001)),
        FieldsAre(2, 0, FloatNear(9.6050, 0.0001)),
        FieldsAre(2, 1, FloatNear(12.4089, 0.0001)),
        FieldsAre(3, 2, FloatNear(7.4041, 0.0001)),
        FieldsAre(4, 0, FloatNear(15.8112, 0.0001)),
        FieldsAre(5, 1, FloatNear(11.6748, 0.0001)),
        FieldsAre(5, 2, FloatNear(11.1788, 0.0001)),
        FieldsAre(5, 3, FloatNear(4.0494, 0.0001)),
        FieldsAre(6, 0, FloatNear(4.5555, 0.0001)),
        FieldsAre(6, 4, FloatNear(8.1384, 0.0001)),
        FieldsAre(7, 0, FloatNear(15.4846, 0.0001)),
        FieldsAre(7, 6, FloatNear(15.9916, 0.0001)),
        FieldsAre(8, 1, FloatNear(14.1860, 0.0001)),
        FieldsAre(8, 2, FloatNear(6.2570, 0.0001)),
        FieldsAre(8, 3, FloatNear(1.3393, 0.0001)),
        FieldsAre(8, 5, FloatNear(1.3845, 0.0001)),
        FieldsAre(9, 0, FloatNear(6.1342, 0.0001)),
        FieldsAre(9, 2, FloatNear(12.5844, 0.0001)),
        FieldsAre(9, 6, FloatNear(2.5074, 0.0001))));
}

TEST(Math, Comparison) {
    EXPECT_TRUE(approximately_equal(95.1, 100.0, 0.05));
    EXPECT_FALSE(essentially_equal(95.1, 100.0, 0.05));

    EXPECT_TRUE(definitely_greater(106.0, 100.0, 0.05));
    EXPECT_FALSE(definitely_greater(95.1, 100.0, 0.05));

    EXPECT_FALSE(definitely_less(105.1, 100.0, 0.05));
    EXPECT_TRUE(definitely_less(94.9, 100.0, 0.05));
}

TEST(Views, SequenceView) {
    // Supporting data
    std::vector<int> vec(5);
    std::iota(vec.begin(), vec.end(), 10);
    std::vector<index_t> indices{1, 3, 4};

    // Only vectors
    SequenceView view(vec, indices);
    EXPECT_THAT(view, ElementsAre(11, 13, 14));
    EXPECT_EQ(view.size(), 3);
    EXPECT_EQ(view[0], 11);
    EXPECT_EQ(view[1], 13);
    EXPECT_EQ(view[2], 14);
    EXPECT_EQ(view.at(0), 11);
    EXPECT_EQ(view.at(1), 13);
    EXPECT_EQ(view.at(2), 14);

    // Vector from SelIndex
    SelIndex sel_indices(indices, 5);
    SequenceView sel_view(vec, sel_indices.indices());
    EXPECT_THAT(sel_view, ElementsAre(11, 13, 14));
    EXPECT_EQ(view.size(), 3);
    EXPECT_EQ(sel_view[0], 11);
    EXPECT_EQ(sel_view[1], 13);
    EXPECT_EQ(sel_view[2], 14);
    EXPECT_EQ(sel_view.at(0), 11);
    EXPECT_EQ(sel_view.at(1), 13);
    EXPECT_EQ(sel_view.at(2), 14);
}
