#ifndef MOLPP_TOOLS_SIMPLEGRAPH_HPP
#define MOLPP_TOOLS_SIMPLEGRAPH_HPP

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

    bool contains_edge(Node const& node1, Node const& node2) const
    {
        auto iter = m_adjacency.find(node1);
        if (iter == m_adjacency.end())
        {
            return false;
        }
        return iter->second.contains(node2);
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

    bool add_edge(Node const& node1, Node const& node2, bool const add_nodes = false)
    {
        // Check nodes existence
        auto adj1_iter = m_adjacency.find(node1);
        if (adj1_iter == m_adjacency.end())
        {
            if (!add_nodes)
            {
                return false;
            }
            adj1_iter = m_adjacency.insert(std::make_pair<Node const&, adjacency_list>(node1, {})).first;
        }

        auto adj2_iter = m_adjacency.find(node2);
        if (adj2_iter == m_adjacency.end())
        {
            if (!add_nodes)
            {
                return false;
            }
            adj2_iter = m_adjacency.insert(std::make_pair<Node const&, adjacency_list>(node2, {})).first;
        }

        adjacency_list& adj1 = adj1_iter->second;
        adjacency_list& adj2 = adj2_iter->second;

        // Insert new edge
        if (!adj1.insert(node2).second)
        {
            return false;
        }
        if (!adj2.insert(node1).second)
        {
            adj1.erase(node2);
            return false;
        }

        return true;
    }

    bool remove_edge(Node const& node1, Node const& node2)
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

        // Remove edge
        bool result = adj1.erase(node2) > 0;
        result &= adj2.erase(node1) > 0;

        return result;
    }

private:
    std::unordered_map<Node, adjacency_list, Hash> m_adjacency;
};

} // namespace mol::internal

#endif // MOLPP_TOOLS_SIMPLEGRAPH_HPP
