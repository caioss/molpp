#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <list>
#include <ranges>
#include <unordered_map>

namespace mol::internal {

template <class Node, class Edge>
class Graph
{
public:
    using node_type = Node;
    using edge_type = Edge;

    Graph()
    {}

    Graph(Graph const &src) = delete;
    Graph &operator=(Graph const &rhs) = delete;

    size_t size() const
    {
        return m_adjacency.size();
    }

    size_t edges_size() const
    {
        return m_edges.size();
    }

    void clear()
    {
        m_adjacency.clear();
        m_edges.clear();
    }

    void clear_edges()
    {
        for (auto &[node, adjacency] : m_adjacency)
        {
            adjacency.clear();
        }

        m_edges.clear();
    }

    bool contains(Node const &node) const
    {
        return m_adjacency.count(node);
    }

    bool add_node(Node const &node)
    {
        return m_adjacency.insert(std::make_pair<Node const &, adjacency_list>(node, {})).second;
    }

    auto const adjacency(Node const &node) const
    {
        return std::views::keys(m_adjacency.at(node));
    }

    auto const nodes() const
    {
        return std::views::keys(m_adjacency);
    }

    auto edges(Node const &node)
    {
        return std::views::values(m_adjacency.at(node)) | std::views::transform([](auto &it){return *it;});
    }

    auto edge_at(Node const &node1, Node const &node2)
    {
        adjacency_list &adj = m_adjacency.at(node1);
        return std::ranges::subrange(adj.find(node2), adj.end()) | std::views::transform([](auto &it){return *(it.second);});
    }

    Edge add_edge(Node const &node1, Node const &node2, Edge const &data)
    {
        auto range = edge_at(node1, node2);
        if (range)
        {
            return *(range.begin());
        }

        // Check nodes existence
        adjacency_list &adj1 = m_adjacency.at(node1);
        adjacency_list &adj2 = m_adjacency.at(node2);

        // Insert new edge
        edge_iterator data_it = m_edges.insert(m_edges.end(), data);
        adj1[node2] = data_it;
        adj2[node1] = data_it;

        return *data_it;
    }

private:
    using edge_list = std::list<Edge>;
    using edge_iterator = typename edge_list::iterator;
    using adjacency_list = std::unordered_map<Node, edge_iterator>;

    // Member variables
    edge_list m_edges;
    std::unordered_map<Node, adjacency_list> m_adjacency;
};

} // namespace mol::internal

#endif // GRAPH_HPP
