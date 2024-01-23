#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <list>
#include <ranges>
#include <optional>
#include <functional>
#include <unordered_map>

namespace mol::internal
{

template<class Node, class EdgeData, class Hash = std::hash<Node>>
class Graph
{
public:
    using node_type = Node;
    using edge_type = EdgeData;
    using edge_result = std::optional<std::reference_wrapper<EdgeData>>;

    Graph() = default;
    ~Graph() = default;
    Graph(Graph&&) = default;
    Graph(Graph const&) = delete;
    Graph& operator=(Graph&&) = default;
    Graph& operator=(Graph const&) = delete;

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
        for (auto& [node, adjacency] : m_adjacency)
        {
            adjacency.clear();
        }

        m_edges.clear();
    }

    bool contains(Node const& node) const
    {
        return m_adjacency.contains(node);
    }

    bool add_node(Node const& node)
    {
        return m_adjacency.insert(std::make_pair<Node const&, adjacency_list>(node, {})).second;
    }

    auto const adjacency(Node const& node) const
    {
        return std::views::keys(m_adjacency.at(node));
    }

    auto const nodes() const
    {
        return std::views::keys(m_adjacency);
    }

    auto edges(Node const& node)
    {
        return std::views::values(m_adjacency.at(node)) | std::views::transform([](auto& it) {
                   return *it;
               });
    }

    edge_result edge_for(Node const& node1, Node const& node2)
    {
        auto adj_iter = m_adjacency.find(node1);
        if (adj_iter == m_adjacency.end())
        {
            return std::nullopt;
        }

        auto edge_iter = adj_iter->second.find(node2);
        if (edge_iter == adj_iter->second.end())
        {
            return std::nullopt;
        }

        return *edge_iter->second;
    }

    edge_result add_edge(Node const& node1, Node const& node2, EdgeData&& data)
    {
        auto edge = edge_for(node1, node2);
        if (edge)
        {
            return edge.value();
        }

        // Check nodes existence
        auto adj1_iter = m_adjacency.find(node1);
        auto adj2_iter = m_adjacency.find(node2);
        if (adj1_iter == m_adjacency.end() || adj2_iter == m_adjacency.end())
        {
            return std::nullopt;
        }
        adjacency_list& adj1 = adj1_iter->second;
        adjacency_list& adj2 = adj2_iter->second;

        // Insert new edge
        edge_iterator data_it = m_edges.insert(m_edges.end(), std::forward<EdgeData>(data));
        adj1[node2] = data_it;
        adj2[node1] = data_it;

        return *data_it;
    }

private:
    using edge_list = std::list<EdgeData>;
    using edge_iterator = typename edge_list::iterator;
    using adjacency_list = std::unordered_map<Node, edge_iterator, Hash>;

    edge_list m_edges;
    std::unordered_map<Node, adjacency_list, Hash> m_adjacency;
};

} // namespace mol::internal

#endif // GRAPH_HPP
