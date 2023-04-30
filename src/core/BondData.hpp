#ifndef BONDGRAPH_HPP
#define BONDGRAPH_HPP

#include <molpp/MolppCore.hpp>
#include "tools/Graph.hpp"
#include <molpp/Bond.hpp>
#include <vector>
#include <memory>
#include <unordered_set>

namespace mol::internal {

class BondData
{
public:
    BondData(size_t const num_atoms);
    BondData() = delete;
    BondData(BondData&&) = default;
    BondData(const BondData &src) = delete;
    BondData &operator=(const BondData &rhs) = delete;

    bool incomplete() const
    {
        return m_incomplete;
    }
    void set_incomplete(bool const incomplete);

    size_t size() const
    {
        return m_graph.edges_size();
    }

    template <class Iterator>
    std::vector<std::shared_ptr<Bond>> bonds(Iterator it, Iterator end);
    template <class Iterator>
    std::vector<index_t> bonded(Iterator it, Iterator end) const;
    std::vector<std::shared_ptr<Bond>> bonds(index_t const index);
    std::shared_ptr<Bond> bond(index_t const atom1, index_t const atom2);
    std::vector<index_t> bonded(index_t const index) const;
    std::shared_ptr<Bond> add_bond(index_t const atom1, index_t const atom2);
    void clear();

private:
    bool m_incomplete;
    Graph<index_t, std::shared_ptr<Bond>> m_graph;
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
std::vector<index_t> BondData::bonded(Iterator it, Iterator end) const
{
    std::unordered_set<index_t> indices;
    while (it != end)
    {
        index_t const index = *(it++);
        auto const range = m_graph.adjacency(index);
        if (range)
        {
            indices.insert(index);
            indices.insert(range.begin(), range.end());
        }
    }
    return std::vector<index_t>(indices.begin(), indices.end());
}

} // namespace mol::internal

#endif // BONDGRAPH_HPP
