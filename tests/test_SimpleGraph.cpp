#include "utils.hpp"
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

// A default constructed graph should have no nodes
TEST_F(SimpleGraphTest, default_constructor)
{
    SimpleGraph<int> default_graph;

    EXPECT_EQ(default_graph.size(), 0);
    EXPECT_THAT(view2vector(default_graph.nodes()), ElementsAre());
}

// SimpleGraph::size should return the number of nodes in the graph
TEST_F(SimpleGraphTest, size)
{
    EXPECT_EQ(graph.size(), num_nodes);
}

// SimpleGraph::clear should remove all nodes from the graph
TEST_F(SimpleGraphTest, clear)
{
    graph.clear();

    EXPECT_EQ(graph.size(), 0);
    EXPECT_THAT(view2vector(graph.nodes()), ElementsAre());
    for (int i = 0; i < num_nodes; i++)
    {
        EXPECT_FALSE(graph.contains(i)) << "Node " << i;
    }
}

// SimpleGraph::contains should return true for nodes that are in the graph
TEST_F(SimpleGraphTest, contains_node)
{
    for (int i = 0; i < num_nodes; i++)
    {
        EXPECT_TRUE(graph.contains(i)) << "Node " << i;
    }
}

// SimpleGraph::contains should return false for nodes that are not in the graph
TEST_F(SimpleGraphTest, contains_invalid_nodes)
{
    EXPECT_FALSE(graph.contains(-1));
    EXPECT_FALSE(graph.contains(num_nodes));
}

// SimpleGraph::contains_edge should return true for edges that are in the graph
TEST_F(SimpleGraphTest, contains_edge)
{
    EXPECT_TRUE(graph.contains_edge(0, 1));
    EXPECT_TRUE(graph.contains_edge(0, 2));
}

// SimpleGraph::contains_edge should return false for edges that are not in the graph
TEST_F(SimpleGraphTest, contains_invalid_edges)
{
    EXPECT_FALSE(graph.contains_edge(0, 0));
    EXPECT_FALSE(graph.contains_edge(0, 3));
    EXPECT_FALSE(graph.contains_edge(1, 2));
}

// SimpleGraph::add_node should add a node to the graph
TEST_F(SimpleGraphTest, add_node)
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

// SimpleGraph::add_node should return false and don't change the graph size when adding an existing node
TEST_F(SimpleGraphTest, add_existing_node)
{
    for (int i = 0; i < num_nodes; i++)
    {
        EXPECT_FALSE(graph.add_node(i)) << "Node " << i;
    }

    EXPECT_EQ(graph.size(), num_nodes);
}

// SimpleGraph::add_node should return false and don't change nodes' adjacency when adding an existing node
TEST_F(SimpleGraphTest, add_existing_node_does_not_change_adjacency)
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

// SimpleGraph::remove_node should remove a node from the graph and return true
TEST_F(SimpleGraphTest, remove_node)
{
    EXPECT_TRUE(graph.remove_node(2));

    EXPECT_FALSE(graph.contains(2));
    EXPECT_THAT(view2vector(graph.nodes()), UnorderedElementsAre(0, 1, 3));
    EXPECT_THAT(view2vector(graph.adjacency(0)), UnorderedElementsAre(1));
    EXPECT_THAT(view2vector(graph.adjacency(1)), UnorderedElementsAre(0));
    EXPECT_THAT(view2vector(graph.adjacency(3)), UnorderedElementsAre());
}

// SimpleGraph::remove_node should return false and don't change the graph when removing an invalid node
TEST_F(SimpleGraphTest, remove_invalid_node)
{
    EXPECT_FALSE(graph.remove_node(4));

    EXPECT_THAT(view2vector(graph.nodes()), UnorderedElementsAre(0, 1, 2, 3));
    EXPECT_THAT(view2vector(graph.adjacency(0)), UnorderedElementsAre(1, 2));
    EXPECT_THAT(view2vector(graph.adjacency(1)), UnorderedElementsAre(0));
    EXPECT_THAT(view2vector(graph.adjacency(2)), UnorderedElementsAre(0));
    EXPECT_THAT(view2vector(graph.adjacency(3)), UnorderedElementsAre());
}

// SimpleGraph::adjacency should return the adjacency list of a node
TEST_F(SimpleGraphTest, adjacency)
{
    EXPECT_THAT(view2vector(graph.adjacency(0)), UnorderedElementsAre(1, 2));
    EXPECT_THAT(view2vector(graph.adjacency(1)), UnorderedElementsAre(0));
    EXPECT_THAT(view2vector(graph.adjacency(2)), UnorderedElementsAre(0));
    EXPECT_THAT(view2vector(graph.adjacency(3)), UnorderedElementsAre());
}

// SimpleGraph::adjacency should throw an exception for invalid nodes
TEST_F(SimpleGraphTest, adjacency_of_invalid_node)
{
    EXPECT_THROW(graph.adjacency(-1), std::out_of_range);
    EXPECT_THROW(graph.adjacency(num_nodes), std::out_of_range);
}

// SimpleGraph::nodes should return the nodes of the graph
TEST_F(SimpleGraphTest, nodes)
{
    EXPECT_THAT(view2vector(graph.nodes()), UnorderedElementsAre(0, 1, 2, 3));
}

// SimpleGraph::clear_edges should remove all edges from the graph
TEST_F(SimpleGraphTest, clear_edges)
{
    graph.clear_edges();

    for (int i = 0; i < num_nodes; i++)
    {
        EXPECT_THAT(view2vector(graph.adjacency(i)), UnorderedElementsAre()) << "Node " << i;
    }
}

// SimpleGraph::add_edge should add an edge to the graph and return true if successful
TEST_F(SimpleGraphTest, add_edge)
{
    EXPECT_TRUE(graph.add_edge(3, 1));
    EXPECT_TRUE(graph.add_edge(3, 2));

    EXPECT_THAT(view2vector(graph.adjacency(0)), UnorderedElementsAre(1, 2));
    EXPECT_THAT(view2vector(graph.adjacency(1)), UnorderedElementsAre(0, 3));
    EXPECT_THAT(view2vector(graph.adjacency(2)), UnorderedElementsAre(0, 3));
    EXPECT_THAT(view2vector(graph.adjacency(3)), UnorderedElementsAre(1, 2));
}

// SimpleGraph::add_edge should return false and don't change the graph when adding an existing edge
TEST_F(SimpleGraphTest, add_existing_edge)
{
    EXPECT_FALSE(graph.add_edge(0, 1));
    EXPECT_FALSE(graph.add_edge(0, 2));

    EXPECT_THAT(view2vector(graph.adjacency(0)), UnorderedElementsAre(1, 2));
    EXPECT_THAT(view2vector(graph.adjacency(1)), UnorderedElementsAre(0));
    EXPECT_THAT(view2vector(graph.adjacency(2)), UnorderedElementsAre(0));
    EXPECT_THAT(view2vector(graph.adjacency(3)), UnorderedElementsAre());
}

// SimpleGraph::add_edge should return false and don't change the graph when adding an edge with invalid nodes
TEST_F(SimpleGraphTest, add_edge_on_invalid_nodes)
{
    EXPECT_FALSE(graph.add_edge(0, -1));
    EXPECT_FALSE(graph.add_edge(-1, 0));
    EXPECT_FALSE(graph.add_edge(-1, -1));
}

// SimpleGraph::add_edge should add both nodes and an edge between them if it's called with add_nodes=true
TEST_F(SimpleGraphTest, add_edge_with_add_nodes)
{
    EXPECT_TRUE(graph.add_edge(4, 5, true));

    EXPECT_THAT(view2vector(graph.nodes()), UnorderedElementsAre(0, 1, 2, 3, 4, 5));
    EXPECT_THAT(view2vector(graph.adjacency(4)), UnorderedElementsAre(5));
}

// SimpleGraph::remove_edge should remove an edge from the graph and return true
TEST_F(SimpleGraphTest, remove_edge)
{
    EXPECT_TRUE(graph.remove_edge(0, 1));

    EXPECT_THAT(view2vector(graph.adjacency(0)), UnorderedElementsAre(2));
    EXPECT_THAT(view2vector(graph.adjacency(1)), UnorderedElementsAre());
    EXPECT_THAT(view2vector(graph.adjacency(2)), UnorderedElementsAre(0));
    EXPECT_THAT(view2vector(graph.adjacency(3)), UnorderedElementsAre());
}

// SimpleGraph::remove_edge should return false and don't change the graph when removing an invalid edge
TEST_F(SimpleGraphTest, remove_invalid_edge)
{
    EXPECT_FALSE(graph.remove_edge(0, 3));

    EXPECT_THAT(view2vector(graph.adjacency(0)), UnorderedElementsAre(1, 2));
    EXPECT_THAT(view2vector(graph.adjacency(1)), UnorderedElementsAre(0));
    EXPECT_THAT(view2vector(graph.adjacency(2)), UnorderedElementsAre(0));
    EXPECT_THAT(view2vector(graph.adjacency(3)), UnorderedElementsAre());
}

// SimpleGraph::remove_edge should return false and don't change the graph when removing an edge with invalid nodes
TEST_F(SimpleGraphTest, remove_edge_on_invalid_nodes)
{
    EXPECT_FALSE(graph.remove_edge(0, -1));
    EXPECT_FALSE(graph.remove_edge(-1, 0));
    EXPECT_FALSE(graph.remove_edge(-1, -1));
}
