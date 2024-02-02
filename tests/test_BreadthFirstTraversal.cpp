#include "tools/algorithms.hpp"
#include <molpp/tools/Graph.hpp>
#include <molpp/tools/SimpleGraph.hpp>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace mol::internal;
using namespace testing;

template<class GraphType>
class BreadthFirstTraversalTest : public testing::Test
{
public:
    BreadthFirstTraversalTest()
    : num_nodes{6}
    , graph()
    , bfs(graph)
    {}

protected:
    void SetUp() override
    {
        for (int i = 0; i < num_nodes; i++)
        {
            ASSERT_TRUE(graph.add_node(i)) << "Node " << i;
        }

        ASSERT_EQ(graph.size(), num_nodes);

        if constexpr (std::same_as<GraphType, Graph<int, int>>)
        {
            ASSERT_TRUE(graph.add_edge(0, 1, 0));
            ASSERT_TRUE(graph.add_edge(0, 2, 0));
            ASSERT_TRUE(graph.add_edge(2, 3, 0));
            ASSERT_TRUE(graph.add_edge(1, 4, 0));
            ASSERT_TRUE(graph.add_edge(4, 5, 0));
            ASSERT_TRUE(graph.add_edge(5, 3, 0));
        }
        else if constexpr (std::same_as<GraphType, SimpleGraph<int>>)
        {
            ASSERT_TRUE(graph.add_edge(0, 1));
            ASSERT_TRUE(graph.add_edge(0, 2));
            ASSERT_TRUE(graph.add_edge(2, 3));
            ASSERT_TRUE(graph.add_edge(1, 4));
            ASSERT_TRUE(graph.add_edge(4, 5));
            ASSERT_TRUE(graph.add_edge(5, 3));
        }
    }

    size_t num_nodes;
    GraphType graph;
    BreadthFirstTraversal<GraphType> bfs;
};

TYPED_TEST_SUITE_P(BreadthFirstTraversalTest);

TYPED_TEST_P(BreadthFirstTraversalTest, InitialState)
{
    EXPECT_THAT(this->bfs.visited(), UnorderedElementsAre());
    EXPECT_THAT(this->bfs.parent_map(), UnorderedElementsAre());
}

TYPED_TEST_P(BreadthFirstTraversalTest, StopAtNode)
{
    EXPECT_TRUE(this->bfs.run(0, [](int const& node) {
        return node == 3;
    }, [](auto){
        return true;
    }));

    // Node 4 may or may not be present (STL-implementation-dependent)
    EXPECT_THAT(this->bfs.visited(), IsSupersetOf({0, 1, 2, 3}));
    EXPECT_THAT(this->bfs.visited(), IsSubsetOf({0, 1, 2, 3, 4}));
    EXPECT_THAT(this->bfs.parent_map(), IsSupersetOf({Pair(3, 2), Pair(2, 0), Pair(1, 0)}));
    EXPECT_THAT(this->bfs.parent_map(), IsSubsetOf({Pair(3, 2), Pair(2, 0), Pair(4, 1), Pair(1, 0)}));
}

TYPED_TEST_P(BreadthFirstTraversalTest, StopAtUnknownNode)
{
    EXPECT_FALSE(this->bfs.run(0, [](int const& node) {
        return node == -1;
    }, [](auto){
        return true;
    }));

    EXPECT_THAT(this->bfs.visited(), UnorderedElementsAre(0, 1, 2, 3, 4, 5));
    // Either edges 4-5 or 5-3 are present (STL-implementation-dependent)
    EXPECT_THAT(this->bfs.parent_map(), IsSupersetOf({Pair(4, 1), Pair(3, 2), Pair(2, 0), Pair(1, 0)}));
    EXPECT_THAT(this->bfs.parent_map(), IsSubsetOf({Pair(5, 3), Pair(5, 4), Pair(4, 1), Pair(3, 2), Pair(2, 0), Pair(1, 0)}));
}

TYPED_TEST_P(BreadthFirstTraversalTest, DontStop)
{
    EXPECT_FALSE(this->bfs.run(0, [](int const& node) {
        return false;
    }, [](auto){
        return true;
    }));

    EXPECT_THAT(this->bfs.visited(), UnorderedElementsAre(0, 1, 2, 3, 4, 5));
    // Either edges 4-5 or 5-3 are present (STL-implementation-dependent)
    EXPECT_THAT(this->bfs.parent_map(), IsSupersetOf({Pair(4, 1), Pair(3, 2), Pair(2, 0), Pair(1, 0)}));
    EXPECT_THAT(this->bfs.parent_map(), IsSubsetOf({Pair(5, 3), Pair(5, 4), Pair(4, 1), Pair(3, 2), Pair(2, 0), Pair(1, 0)}));
}

TYPED_TEST_P(BreadthFirstTraversalTest, StopAtNodeWithMask)
{
    EXPECT_TRUE(this->bfs.run(0, [](int const& node) {
        return node == 3;
    }, [](int const& node){
        return node != 4;
    }));

    EXPECT_THAT(this->bfs.visited(), UnorderedElementsAre(0, 1, 2, 3));
    EXPECT_THAT(this->bfs.parent_map(), UnorderedElementsAre(Pair(3, 2), Pair(2, 0), Pair(1, 0)));
}

TYPED_TEST_P(BreadthFirstTraversalTest, MaskAllNodes)
{
    EXPECT_FALSE(this->bfs.run(0, [](int const& node) {
        return false;
    }, [](int const& node){
        return false;
    }));

    EXPECT_THAT(this->bfs.visited(), UnorderedElementsAre());
    EXPECT_THAT(this->bfs.parent_map(), UnorderedElementsAre());
}

REGISTER_TYPED_TEST_SUITE_P(BreadthFirstTraversalTest, InitialState, StopAtNode, StopAtUnknownNode, DontStop, StopAtNodeWithMask, MaskAllNodes);

using GraphTypes = testing::Types<Graph<int, int>, SimpleGraph<int>>;
INSTANTIATE_TYPED_TEST_SUITE_P(BreadthFirstTraversalAllGraphs, BreadthFirstTraversalTest, GraphTypes);
