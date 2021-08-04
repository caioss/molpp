#include "BondGraph.hpp"

using namespace mol;
using namespace mol::internal;

BondGraph::BondGraph(size_t const num_atoms)
: m_incomplete { true }
{
    for (size_t i = 0; i < num_atoms; i++)
    {
        m_graph.add_node(i);
    }
}

void BondGraph::set_incomplete(bool const incomplete)
{
    m_incomplete = incomplete;
}

std::vector<std::shared_ptr<Bond>> BondGraph::bonds(size_t const index)
{
    auto range = m_graph.edges(index);
    return {range.begin(), range.end()};
}

std::shared_ptr<Bond> BondGraph::bond(size_t const atom1, size_t const atom2)
{
    auto range = m_graph.edge_at(atom1, atom2);
    if (!range.is_valid())
    {
        return nullptr;
    }
    return *(range.begin());
}

std::vector<size_t> BondGraph::bonded(size_t const index) const
{
    auto const range = m_graph.adjacency(index);
    return {range.begin(), range.end()};
}

std::shared_ptr<Bond> BondGraph::add_bond(size_t const atom1, size_t const atom2)
{
    std::shared_ptr<Bond> data = std::make_shared<Bond>(atom1, atom2);
    return *m_graph.add_edge(atom1, atom2, data);
}
