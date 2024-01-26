#include "auxiliary.hpp"
#include <molpp/tools/SimpleGraph.hpp>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <vector>
#include <numeric>

using namespace mol::internal;
using namespace testing;

class SimpleGraphTest : public testing::Test
{
public:
    SimpleGraphTest()
    : num_nodes{4}
    {}

protected:
    void SetUp() override
    {
        for (int i = 0; i < num_nodes; i++)
        {
            ASSERT_TRUE(graph.add_node(i)) << "Node " << i;
        }
        ASSERT_TRUE(graph.add_edge(0, 1));
        ASSERT_TRUE(graph.add_edge(0, 2));

        ASSERT_EQ(graph.size(), num_nodes);
    }

    size_t num_nodes;
    SimpleGraph<int> graph;
};

TEST_F(SimpleGraphTest, DefaultConstructor)
{
    SimpleGraph<int> default_graph;

    EXPECT_EQ(default_graph.size(), 0);
    EXPECT_THAT(view2vector(default_graph.nodes()), ElementsAre());
}

TEST_F(SimpleGraphTest, Size)
{
    EXPECT_EQ(graph.size(), num_nodes);
}

TEST_F(SimpleGraphTest, Clear)
{
    graph.clear();

    EXPECT_EQ(graph.size(), 0);
    EXPECT_THAT(view2vector(graph.nodes()), ElementsAre());
    for (int i = 0; i < num_nodes; i++)
    {
        EXPECT_FALSE(graph.contains(i)) << "Node " << i;
    }
}

TEST_F(SimpleGraphTest, ContainsNode)
{
    for (int i = 0; i < num_nodes; i++)
    {
        EXPECT_TRUE(graph.contains(i)) << "Node " << i;
    }
}

TEST_F(SimpleGraphTest, ContainsInvalidNodes)
{
    EXPECT_FALSE(graph.contains(-1));
    EXPECT_FALSE(graph.contains(num_nodes));
}

TEST_F(SimpleGraphTest, AddNode)
{
    size_t const extra_nodes = 3;
    for (int i = num_nodes; i < num_nodes + extra_nodes; i++)
    {
        EXPECT_TRUE(graph.add_node(i)) << "Node " << i;
        EXPECT_TRUE(graph.contains(i)) << "Node " << i;
        EXPECT_THAT(view2vector(graph.adjacency(i)), UnorderedElementsAre()) << "Node " << i;
    }
    EXPECT_EQ(graph.size(), num_nodes + extra_nodes);
    EXPECT_THAT(view2vector(graph.nodes()), UnorderedElementsAre(0, 1, 2, 3, 4, 5, 6));
}

TEST_F(SimpleGraphTest, AddExistingNode)
{
    for (int i = 0; i < num_nodes; i++)
    {
        EXPECT_FALSE(graph.add_node(i)) << "Node " << i;
    }

    EXPECT_EQ(graph.size(), num_nodes);
}

TEST_F(SimpleGraphTest, AddExistingNodeDoesNotChangeAdjacency)
{
    for (int i = 0; i < num_nodes; i++)
    {
        EXPECT_FALSE(graph.add_node(i)) << "Node " << i;
    }

    EXPECT_THAT(view2vector(graph.adjacency(0)), UnorderedElementsAre(1, 2));
    EXPECT_THAT(view2vector(graph.adjacency(1)), UnorderedElementsAre(0));
    EXPECT_THAT(view2vector(graph.adjacency(2)), UnorderedElementsAre(0));
    EXPECT_THAT(view2vector(graph.adjacency(3)), UnorderedElementsAre());
}

TEST_F(SimpleGraphTest, Adjacency)
{
    EXPECT_THAT(view2vector(graph.adjacency(0)), UnorderedElementsAre(1, 2));
    EXPECT_THAT(view2vector(graph.adjacency(1)), UnorderedElementsAre(0));
    EXPECT_THAT(view2vector(graph.adjacency(2)), UnorderedElementsAre(0));
    EXPECT_THAT(view2vector(graph.adjacency(3)), UnorderedElementsAre());
}

TEST_F(SimpleGraphTest, AdjacencyOfInvalidNode)
{
    EXPECT_THROW(graph.adjacency(-1), std::out_of_range);
    EXPECT_THROW(graph.adjacency(num_nodes), std::out_of_range);
}

TEST_F(SimpleGraphTest, Nodes)
{
    EXPECT_THAT(view2vector(graph.nodes()), UnorderedElementsAre(0, 1, 2, 3));
}

TEST_F(SimpleGraphTest, ClearEdges)
{
    graph.clear_edges();

    for (int i = 0; i < num_nodes; i++)
    {
        EXPECT_THAT(view2vector(graph.adjacency(i)), UnorderedElementsAre()) << "Node " << i;
    }
}

TEST_F(SimpleGraphTest, AddEdge)
{
    EXPECT_TRUE(graph.add_edge(3, 1));
    EXPECT_TRUE(graph.add_edge(3, 2));

    EXPECT_THAT(view2vector(graph.adjacency(0)), UnorderedElementsAre(1, 2));
    EXPECT_THAT(view2vector(graph.adjacency(1)), UnorderedElementsAre(0, 3));
    EXPECT_THAT(view2vector(graph.adjacency(2)), UnorderedElementsAre(0, 3));
    EXPECT_THAT(view2vector(graph.adjacency(3)), UnorderedElementsAre(1, 2));
}

TEST_F(SimpleGraphTest, AddExistingEdge)
{
    EXPECT_TRUE(graph.add_edge(0, 1));
    EXPECT_TRUE(graph.add_edge(0, 2));

    EXPECT_THAT(view2vector(graph.adjacency(0)), UnorderedElementsAre(1, 2));
    EXPECT_THAT(view2vector(graph.adjacency(1)), UnorderedElementsAre(0));
    EXPECT_THAT(view2vector(graph.adjacency(2)), UnorderedElementsAre(0));
    EXPECT_THAT(view2vector(graph.adjacency(3)), UnorderedElementsAre());
}

TEST_F(SimpleGraphTest, AddEdgeOnInvalidNodes)
{
    EXPECT_FALSE(graph.add_edge(0, -1));
    EXPECT_FALSE(graph.add_edge(-1, 0));
    EXPECT_FALSE(graph.add_edge(-1, -1));
}
