#ifndef BONDGRAPH_HPP
#define BONDGRAPH_HPP

#include "Bond.hpp"
#include <vector>
#include <memory>
#include <unordered_set>
#include <unordered_map>

namespace mol::internal {

class BondGraph
{
public:
    BondGraph() = delete;
    BondGraph(size_t const num_atoms);
    BondGraph(const BondGraph &src) = delete;
    BondGraph &operator=(const BondGraph &rhs) = delete;

    bool incomplete() const { return m_incomplete; }
    void set_incomplete(bool const incomplete);

    size_t size() const { return m_graph.size(); }

    template <class Iterator>
    std::vector<std::shared_ptr<Bond>> bonds(Iterator it, Iterator end);
    template <class Iterator>
    std::vector<size_t> bonded(Iterator it, Iterator end) const;
    std::vector<std::shared_ptr<Bond>> bonds(size_t const index);
    std::shared_ptr<Bond> bond(size_t const atom1, size_t const atom2);
    std::vector<size_t> bonded(size_t const index) const;
    std::shared_ptr<Bond> add_bond(size_t const atom1, size_t const atom2);

private:
    bool m_incomplete;
    std::vector<std::unordered_map<size_t, std::shared_ptr<Bond>>> m_graph;
};

template <class Iterator>
std::vector<std::shared_ptr<Bond>> BondGraph::bonds(Iterator it, Iterator end)
{
    std::unordered_set<std::shared_ptr<Bond>> indices;
    while (it != end)
    {
        for (auto const &item : m_graph[*(it++)])
        {
            indices.insert(item.second);
        }
    }
    return std::vector<std::shared_ptr<Bond>>(indices.begin(), indices.end());
}

template <class Iterator>
std::vector<size_t> BondGraph::bonded(Iterator it, Iterator end) const
{
    std::unordered_set<size_t> indices;
    while (it != end)
    {
        size_t const index = *(it++);
        auto const &adjacency = m_graph[index];
        if (adjacency.size())
        {
            indices.insert(index);
        }
        for (auto const &item : adjacency)
        {
            indices.insert(item.first);
        }
    }
    return std::vector<size_t>(indices.begin(), indices.end());
}

} // namespace mol::internal

#endif // BONDGRAPH_HPP
