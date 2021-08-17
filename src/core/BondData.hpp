#ifndef BONDGRAPH_HPP
#define BONDGRAPH_HPP

#include "tools/Graph.hpp"
#include <molpp/Bond.hpp>
#include <vector>
#include <memory>
#include <unordered_set>

namespace mol::internal {

class BondData
{
public:
    BondData() = delete;
    BondData(size_t const num_atoms);
    BondData(const BondData &src) = delete;
    BondData &operator=(const BondData &rhs) = delete;

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
    Graph<size_t, std::shared_ptr<Bond>> m_graph;
};

template <class Iterator>
std::vector<std::shared_ptr<Bond>> BondData::bonds(Iterator it, Iterator end)
{
    std::unordered_set<std::shared_ptr<Bond>> indices;
    while (it != end)
    {
        auto range = m_graph.edges(*(it++));
        indices.insert(range.begin(), range.end());
    }
    return std::vector<std::shared_ptr<Bond>>(indices.begin(), indices.end());
}

template <class Iterator>
std::vector<size_t> BondData::bonded(Iterator it, Iterator end) const
{
    std::unordered_set<size_t> indices;
    while (it != end)
    {
        size_t const index = *(it++);
        auto const range = m_graph.adjacency(index);
        if (range.is_valid())
        {
            indices.insert(index);
            indices.insert(range.begin(), range.end());
        }
    }
    return std::vector<size_t>(indices.begin(), indices.end());
}

} // namespace mol::internal

#endif // BONDGRAPH_HPP
