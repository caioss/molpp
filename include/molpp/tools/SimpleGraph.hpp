#ifndef SIMPLEGRAPH_HPP
#define SIMPLEGRAPH_HPP

#include <ranges>
#include <unordered_map>
#include <unordered_set>

namespace mol::internal
{

template<class Node, class Hash = std::hash<Node>>
class SimpleGraph
{
public:
    using node_type = Node;
    using adjacency_list = std::unordered_set<Node, Hash>;

    SimpleGraph() = default;
    ~SimpleGraph() = default;
    SimpleGraph(SimpleGraph&&) = default;
    SimpleGraph(SimpleGraph const&) = delete;
    SimpleGraph& operator=(SimpleGraph&&) = default;
    SimpleGraph& operator=(SimpleGraph const&) = delete;

    size_t size() const
    {
        return m_adjacency.size();
    }

    void clear()
    {
        m_adjacency.clear();
    }

    void clear_edges()
    {
        for (auto& [node, adjacency] : m_adjacency)
        {
            adjacency.clear();
        }
    }

    bool contains(Node const& node) const
    {
        return m_adjacency.contains(node);
    }

    bool add_node(Node const& node)
    {
        return m_adjacency.insert(std::make_pair<Node const&, adjacency_list>(node, {})).second;
    }

    bool remove_node(Node const& node)
    {
        // Check for node existence
        auto iter = m_adjacency.find(node);
        if (iter == m_adjacency.end())
        {
            return false;
        }
        adjacency_list& adjacency = iter->second;

        // Remove edges
        for (Node const& adj_node : adjacency)
        {
            m_adjacency[adj_node].erase(node);
        }

        // Remove the node itself
        m_adjacency.erase(iter);

        return true;
    }

    adjacency_list const& adjacency(Node const& node) const
    {
        return m_adjacency.at(node);
    }

    auto const nodes() const
    {
        return std::views::keys(m_adjacency);
    }

    bool add_edge(Node const& node1, Node const& node2)
    {
        // Check nodes existence
        auto adj1_iter = m_adjacency.find(node1);
        auto adj2_iter = m_adjacency.find(node2);
        if (adj1_iter == m_adjacency.end() || adj2_iter == m_adjacency.end())
        {
            return false;
        }
        adjacency_list& adj1 = adj1_iter->second;
        adjacency_list& adj2 = adj2_iter->second;

        // Insert new edge
        adj1.insert(node2);
        adj2.insert(node1);

        return true;
    }

private:
    std::unordered_map<Node, adjacency_list, Hash> m_adjacency;
};

} // namespace mol::internal

#endif // SIMPLEGRAPH_HPP
