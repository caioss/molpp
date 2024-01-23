#include <molpp/internal/BondData.hpp>

using namespace mol;
using namespace mol::internal;

BondData::BondData(size_t const num_atoms)
: m_incomplete { true }
{
    for (index_t i = 0; i < num_atoms; i++)
    {
        m_graph.add_node(i);
    }
}

void BondData::set_incomplete(bool const incomplete)
{
    m_incomplete = incomplete;
}

std::vector<std::shared_ptr<Bond>> BondData::bonds(index_t const index)
{
    auto range = m_graph.edges(index);
    return {range.begin(), range.end()};
}

std::shared_ptr<Bond> BondData::bond(index_t const atom1, index_t const atom2)
{
    auto edge = m_graph.edge_for(atom1, atom2);
    if (!edge)
    {
        return nullptr;
    }
    return edge.value();
}

std::vector<index_t> BondData::bonded(index_t const index) const
{
    auto const range = m_graph.adjacency(index);
    if (!range)
    {
        return {};
    }

    std::vector<index_t> indices{range.begin(), range.end()};
    indices.push_back(index);
    return indices;
}

std::shared_ptr<Bond> BondData::add_bond(index_t const atom1, index_t const atom2)
{
    if (atom1 == atom2)
    {
        return nullptr;
    }

    auto result = m_graph.add_edge(atom1, atom2, std::make_shared<Bond>(atom1, atom2));
    return result ? result.value().get() : nullptr;
}

void BondData::clear()
{
    m_graph.clear_edges();
}
