#ifndef GRAPH_HPP
#define GRAPH_HPP

#include "tools/iterators.hpp"
#include <list>
#include <unordered_map>

namespace mol::internal {

template <class Node, class Edge>
class Graph
{
private:
    using edge_list = std::list<Edge>;
    using edge_iterator = typename edge_list::iterator;
    using adjacency_list = std::unordered_map<Node, edge_iterator>;

    template <class It, class Type>
    class BaseNodeIt;

    template <class It, class Type>
    class BaseEdgeIt;

    // Member variables
    edge_list m_edges;
    std::unordered_map<Node, adjacency_list> m_adjacency;

public:
    using node_type = Node;
    using edge_type = Edge;
    using NodeIt = BaseNodeIt<typename adjacency_list::iterator, Node>;
    using ConstNodeIt = BaseNodeIt<typename adjacency_list::const_iterator, const Node>;
    using EdgeIt = BaseEdgeIt<typename adjacency_list::iterator, Edge>;
    using ConstEdgeIt = BaseEdgeIt<typename adjacency_list::const_iterator, const Edge>;

    Graph()
    {}

    Graph(const Graph &src) = delete;
    Graph &operator=(const Graph &rhs) = delete;

    size_t size() const {
        return m_adjacency.size();
    }

    bool contains(Node const &node) const
    {
        return m_adjacency.count(node);
    }

    bool add_node(Node const &node)
    {
        return m_adjacency.insert(std::make_pair<Node const &, adjacency_list>(node, {})).second;
    }

    Range<ConstNodeIt> adjacency(Node const &node) const
    {
        adjacency_list const &adj = m_adjacency.at(node);
        return Range(ConstNodeIt(adj.cbegin()), ConstNodeIt(adj.cend()));
    }

    Range<EdgeIt> edges(Node const &node)
    {
        adjacency_list &adj = m_adjacency.at(node);
        return Range(EdgeIt(adj.begin()), EdgeIt(adj.end()));
    }

    Range<EdgeIt> edge_at(Node const &node1, Node const &node2)
    {
        adjacency_list &adj = m_adjacency.at(node1);
        return Range(EdgeIt(adj.find(node2)), EdgeIt(adj.end()));
    }

    EdgeIt add_edge(Node const &node1, Node const &node2, Edge const &data)
    {
        Range<EdgeIt> edge_it = edge_at(node1, node2);
        if (edge_it.is_valid())
        {
            return edge_it.begin();
        }

        // Check nodes existence
        adjacency_list &adj1 = m_adjacency.at(node1);
        adjacency_list &adj2 = m_adjacency.at(node2);

        // Insert new edge
        edge_iterator data_it = m_edges.insert(m_edges.end(), data);
        adj2[node1] = data_it;
        auto iterator = adj1.emplace(node2, data_it).first;

        return EdgeIt(iterator);
    }

private:
    // Iterators
    template <class It, class Type>
    class BaseNodeIt : public IteratorWrapper<It>
    {
        using base = IteratorWrapper<It>;
        using base::m_iterator;

    public:
        using iterator_category = typename base::iterator_category;
        using difference_type = typename base::difference_type;
        using value_type = Type;
        using pointer = value_type *;
        using reference = value_type &;

        using base::IteratorWrapper;

        reference operator*() const
        {
            return (*m_iterator).first;
        }
    };

    template <class It, class Type>
    class BaseEdgeIt : public IteratorWrapper<It>
    {
        using base = IteratorWrapper<It>;
        using base::m_iterator;

    public:
        using iterator_category = typename base::iterator_category;
        using difference_type = typename base::difference_type;
        using value_type = Type;
        using pointer = value_type *;
        using reference = value_type &;

        using base::IteratorWrapper;

        reference operator*() const
        {
            return *((*m_iterator).second);
        }
    };
};

} // namespace mol::internal

#endif // GRAPH_HPP
