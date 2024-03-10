#include "auxiliary.hpp"
#include <molpp/tools/Graph.hpp>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <vector>
#include <numeric>

using namespace mol::internal;
using namespace testing;

class GraphTest : public testing::Test
{
public:
    GraphTest()
    : num_nodes{4}
    , num_edges{2}
    {}

protected:
    void SetUp() override
    {
        for (size_t i = 0; i < num_nodes; i++)
        {
            ASSERT_TRUE(graph.add_node(i)) << "Node " << i;
        }
        ASSERT_EQ(graph.add_edge(0, 1, -1), -1);
        ASSERT_EQ(graph.add_edge(0, 2, -2), -2);

        ASSERT_EQ(graph.size(), num_nodes);
        ASSERT_EQ(graph.edges_size(), num_edges);
    }

    size_t num_nodes;
    size_t num_edges;
    Graph<int, int> graph;
};

TEST_F(GraphTest, DefaultConstructor)
{
    Graph<int, int> default_graph;

    EXPECT_EQ(default_graph.size(), 0);
    EXPECT_EQ(default_graph.edges_size(), 0);
    EXPECT_THAT(view2vector(default_graph.nodes()), ElementsAre());
}

TEST_F(GraphTest, Clear)
{
    graph.clear();

    EXPECT_EQ(graph.size(), 0);
    EXPECT_EQ(graph.edges_size(), 0);
    EXPECT_THAT(view2vector(graph.nodes()), ElementsAre());
    for (size_t i = 0; i < num_nodes; i++)
    {
        EXPECT_FALSE(graph.contains(i)) << "Node " << i;
    }
}

TEST_F(GraphTest, ContainsNode)
{
    for (size_t i = 0; i < num_nodes; i++)
    {
        EXPECT_TRUE(graph.contains(i)) << "Node " << i;
    }
}

TEST_F(GraphTest, ContainsInvalidNodes)
{
    EXPECT_FALSE(graph.contains(-1));
    EXPECT_FALSE(graph.contains(num_nodes));
}

TEST_F(GraphTest, AddNode)
{
    size_t const extra_nodes = 3;
    for (size_t i = num_nodes; i < num_nodes + extra_nodes; i++)
    {
        EXPECT_TRUE(graph.add_node(i)) << "Node " << i;
        EXPECT_TRUE(graph.contains(i)) << "Node " << i;
        EXPECT_THAT(view2vector(graph.adjacency(i)), UnorderedElementsAre()) << "Node " << i;
    }
    EXPECT_EQ(graph.size(), num_nodes + extra_nodes);
    EXPECT_THAT(view2vector(graph.nodes()), UnorderedElementsAre(0, 1, 2, 3, 4, 5, 6));
}

TEST_F(GraphTest, AddExistingNode)
{
    for (size_t i = 0; i < num_nodes; i++)
    {
        EXPECT_FALSE(graph.add_node(i)) << "Node " << i;
    }

    EXPECT_EQ(graph.size(), num_nodes);
}

TEST_F(GraphTest, AddExistingNodeDoesNotChangeEdgeData)
{
    for (size_t i = 0; i < num_nodes; i++)
    {
        EXPECT_FALSE(graph.add_node(i)) << "Node " << i;
    }

    EXPECT_THAT(view2vector(graph.edges(0)), UnorderedElementsAre(-1, -2));
    EXPECT_THAT(view2vector(graph.edges(1)), UnorderedElementsAre(-1));
    EXPECT_THAT(view2vector(graph.edges(2)), UnorderedElementsAre(-2));
}

TEST_F(GraphTest, AddExistingNodeDoesNotChangeAdjacency)
{
    for (size_t i = 0; i < num_nodes; i++)
    {
        EXPECT_FALSE(graph.add_node(i)) << "Node " << i;
    }

    EXPECT_THAT(view2vector(graph.adjacency(0)), UnorderedElementsAre(1, 2));
    EXPECT_THAT(view2vector(graph.adjacency(1)), UnorderedElementsAre(0));
    EXPECT_THAT(view2vector(graph.adjacency(2)), UnorderedElementsAre(0));
    EXPECT_THAT(view2vector(graph.adjacency(3)), UnorderedElementsAre());
}

TEST_F(GraphTest, Adjacency)
{
    EXPECT_THAT(view2vector(graph.adjacency(0)), UnorderedElementsAre(1, 2));
    EXPECT_THAT(view2vector(graph.adjacency(1)), UnorderedElementsAre(0));
    EXPECT_THAT(view2vector(graph.adjacency(2)), UnorderedElementsAre(0));
    EXPECT_THAT(view2vector(graph.adjacency(3)), UnorderedElementsAre());
}

TEST_F(GraphTest, AdjacencyOfInvalidNode)
{
    EXPECT_THROW(graph.adjacency(-1), std::out_of_range);
    EXPECT_THROW(graph.adjacency(num_nodes), std::out_of_range);
}

TEST_F(GraphTest, Nodes)
{
    EXPECT_THAT(view2vector(graph.nodes()), UnorderedElementsAre(0, 1, 2, 3));
}

TEST_F(GraphTest, Edges)
{
    EXPECT_THAT(view2vector(graph.edges(0)), UnorderedElementsAre(-1, -2));
    EXPECT_THAT(view2vector(graph.edges(1)), UnorderedElementsAre(-1));
    EXPECT_THAT(view2vector(graph.edges(2)), UnorderedElementsAre(-2));
    EXPECT_THAT(view2vector(graph.edges(3)), UnorderedElementsAre());
}

TEST_F(GraphTest, EdgesOfInvalidNode)
{
    EXPECT_THROW(graph.edges(-1), std::out_of_range);
    EXPECT_THROW(graph.edges(num_nodes), std::out_of_range);
}

TEST_F(GraphTest, EdgeFor)
{
    EXPECT_EQ(graph.edge_for(0, 1).value(), -1);
    EXPECT_EQ(graph.edge_for(0, 2).value(), -2);
    EXPECT_FALSE(graph.edge_for(0, 3));
}

TEST_F(GraphTest, EdgeForWithInavlidNode)
{
    EXPECT_FALSE(graph.edge_for(-1, 0));
    EXPECT_FALSE(graph.edge_for(0, -1));
}

TEST_F(GraphTest, ClearEdges)
{
    graph.clear_edges();

    EXPECT_FALSE(graph.edge_for(0, 1));
    EXPECT_FALSE(graph.edge_for(0, 2));
}

TEST_F(GraphTest, AddEdges)
{
    EXPECT_EQ(graph.add_edge(3, 1, -1), -1);
    EXPECT_EQ(graph.add_edge(3, 2, -2), -2);

    EXPECT_EQ(graph.edge_for(3, 1).value(), -1);
    EXPECT_EQ(graph.edge_for(3, 2).value(), -2);
}

TEST_F(GraphTest, AddExistingEdge)
{
    graph.clear_edges();
    int& first_edge = graph.add_edge(0, 1, -1).value();
    int& second_edge = graph.add_edge(0, 1, -1).value();

    EXPECT_EQ(&first_edge, &second_edge);
}

TEST_F(GraphTest, AddEdgeOnInvalidNodes)
{
    EXPECT_FALSE(graph.add_edge(0, -1, -1));
    EXPECT_FALSE(graph.add_edge(-1, 0, -1));
    EXPECT_FALSE(graph.add_edge(-1, -1, -1));
}
