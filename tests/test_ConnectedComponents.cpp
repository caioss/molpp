#include "tools/algorithms.hpp"
#include <molpp/tools/Graph.hpp>
#include <molpp/tools/SimpleGraph.hpp>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace mol::internal;
using namespace testing;

template<class GraphType>
class ConnectedComponentsTest : public testing::Test
{
public:
    ConnectedComponentsTest()
    : num_nodes{7}
    , graph()
    , components(graph)
    {}

protected:
    void SetUp() override
    {
        for (size_t i = 0; i < num_nodes; i++)
        {
            ASSERT_TRUE(graph.add_node(i)) << "Node " << i;
        }

        ASSERT_EQ(graph.size(), num_nodes);

        if constexpr (std::same_as<GraphType, Graph<int, int>>)
        {
            ASSERT_TRUE(graph.add_edge(0, 1, 0));
            ASSERT_TRUE(graph.add_edge(0, 2, 0));
            ASSERT_TRUE(graph.add_edge(2, 3, 0));
            ASSERT_TRUE(graph.add_edge(4, 5, 0));
        }
        else if constexpr (std::same_as<GraphType, SimpleGraph<int>>)
        {
            ASSERT_TRUE(graph.add_edge(0, 1));
            ASSERT_TRUE(graph.add_edge(0, 2));
            ASSERT_TRUE(graph.add_edge(2, 3));
            ASSERT_TRUE(graph.add_edge(4, 5));
        }
    }

    size_t num_nodes;
    GraphType graph;
    ConnectedComponents<GraphType> components;
};

TYPED_TEST_SUITE_P(ConnectedComponentsTest);

TYPED_TEST_P(ConnectedComponentsTest, InitialState)
{
    EXPECT_THAT(this->components.components(), UnorderedElementsAre());
}

TYPED_TEST_P(ConnectedComponentsTest, WithoutEdges)
{
    this->graph.clear_edges();

    EXPECT_EQ(this->components.run([](auto){return true;}), 7);

    EXPECT_THAT(this->components.components(), UnorderedElementsAre(ElementsAre(0), ElementsAre(1), ElementsAre(2), ElementsAre(3), ElementsAre(4), ElementsAre(5), ElementsAre(6)));
}

TYPED_TEST_P(ConnectedComponentsTest, WithoutFilter)
{
    EXPECT_EQ(this->components.run([](auto){return true;}), 3);

    EXPECT_THAT(this->components.components(), UnorderedElementsAre(UnorderedElementsAre(0, 1, 2, 3), UnorderedElementsAre(4, 5), ElementsAre(6)));
}

TYPED_TEST_P(ConnectedComponentsTest, WithFilter)
{
    EXPECT_EQ(this->components.run([](auto const &node){return node != 2;}), 4);

    EXPECT_THAT(this->components.components(), UnorderedElementsAre(UnorderedElementsAre(0, 1), ElementsAre(3), UnorderedElementsAre(4, 5), ElementsAre(6)));
}

TYPED_TEST_P(ConnectedComponentsTest, FilterAll)
{
    EXPECT_EQ(this->components.run([](auto){return false;}), 0);

    EXPECT_THAT(this->components.components(), ElementsAre());
}

REGISTER_TYPED_TEST_SUITE_P(ConnectedComponentsTest, InitialState, WithoutEdges, WithoutFilter, WithFilter, FilterAll);

using GraphTypes = testing::Types<Graph<int, int>, SimpleGraph<int>>;
INSTANTIATE_TYPED_TEST_SUITE_P(ConnectedComponentsAllGraphs, ConnectedComponentsTest, GraphTypes);
