#include "tools/iterators.hpp"
#include "tools/Graph.hpp"
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <vector>

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

TEST(Iterators, Range) {
    using iter_type = typename std::vector<int>::iterator;
    using RangeInt = mol::internal::Range<iter_type>;
    std::vector<int> values(1);

    RangeInt range(values.begin(), values.end());
    EXPECT_EQ(range.begin(), values.begin());
    EXPECT_EQ(range.end(), values.end());
    EXPECT_TRUE(range.is_valid());

    RangeInt invalid(values.end(), values.end());
    EXPECT_EQ(invalid.begin(), values.end());
    EXPECT_EQ(invalid.end(), values.end());
    EXPECT_FALSE(invalid.is_valid());
}

TEST(Graph, Graph) {
    Graph<int, int> graph;

    EXPECT_EQ(graph.size(), 0);
    EXPECT_TRUE(graph.add_node(0));
    EXPECT_FALSE(graph.add_node(0));
    EXPECT_TRUE(graph.add_node(1));
    EXPECT_TRUE(graph.add_node(2));
    EXPECT_TRUE(graph.add_node(3));
    EXPECT_EQ(graph.size(), 4);
    EXPECT_TRUE(graph.contains(0));
    EXPECT_FALSE(graph.contains(-1));

    EXPECT_EQ(*graph.add_edge(0, 1, -1), -1);
    EXPECT_EQ(*graph.add_edge(0, 2, -2), -2);

    EXPECT_THAT(graph.adjacency(0), UnorderedElementsAre(1, 2));
    EXPECT_THAT(graph.adjacency(1), UnorderedElementsAre(0));
    EXPECT_THAT(graph.adjacency(2), UnorderedElementsAre(0));
    EXPECT_THAT(graph.adjacency(3), UnorderedElementsAre());

    EXPECT_THAT(graph.edges(0), UnorderedElementsAre(-1, -2));
    EXPECT_THAT(graph.edges(1), UnorderedElementsAre(-1));
    EXPECT_THAT(graph.edges(2), UnorderedElementsAre(-2));
    EXPECT_THAT(graph.edges(3), UnorderedElementsAre());

    EXPECT_THAT(graph.edge_at(1, 0), UnorderedElementsAre(-1));
    EXPECT_THAT(graph.edge_at(2, 1), UnorderedElementsAre());
}
