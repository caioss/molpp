#include "BondGraph.hpp"

using namespace mol;
using namespace mol::internal;

BondGraph::BondGraph(size_t const num_atoms)
: m_incomplete { true },
  m_graph(num_atoms)
{
}

void BondGraph::set_incomplete(bool const incomplete)
{
    m_incomplete = incomplete;
}

std::vector<std::shared_ptr<Bond>> BondGraph::bonds(size_t const index)
{
    auto &adjacency = m_graph[index];
    std::vector<std::shared_ptr<Bond>> bonded;
    bonded.reserve(adjacency.size());
    for (auto const &item : adjacency)
    {
        bonded.push_back(item.second);
    }
    return bonded;
}

std::shared_ptr<Bond> BondGraph::bond(size_t const atom1, size_t const atom2)
{
    if (!m_graph[atom1].count(atom2))
    {
        return nullptr;
    }
    // Note: operator[] may insert items.
    // We're good because the presence of atom2 was already tested.
    return m_graph[atom1][atom2];
}

std::vector<size_t> BondGraph::bonded(size_t const index) const
{
    std::vector<size_t> indices;
    indices.reserve(m_graph[index].size());
    for (auto const &item : m_graph[index])
    {
        indices.push_back(item.first);
    }
    return indices;
}

std::shared_ptr<Bond> BondGraph::add_bond(size_t const atom1, size_t const atom2)
{
    auto &atom1_adj = m_graph[atom1];
    auto &atom2_adj = m_graph[atom2];
    if (atom1_adj.count(atom2) || atom2_adj.count(atom1))
    {
        return nullptr;
    }

    std::shared_ptr<Bond> data = std::make_shared<Bond>(atom1, atom2);
    if (!atom1_adj.insert(std::make_pair(atom2, data)).second)
    {
        return nullptr;
    }
    if (!atom2_adj.insert(std::make_pair(atom1, data)).second)
    {
        atom1_adj.erase(atom2);
        return nullptr;
    }

    return data;
}
