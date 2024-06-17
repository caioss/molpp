#include <molpp/internal/Topology.hpp>

namespace mol::internal
{

bool operator==(Topology::MolecularEntityId const& lhs, Topology::MolecularEntityId const& rhs)
{
    return lhs.category == rhs.category && lhs.index == rhs.index;
}

bool Topology::add_link(MolecularEntityId const entity1, MolecularEntityId const entity2)
{
    return m_topology.add_edge(entity1, entity2, true);
}

bool Topology::remove_link(MolecularEntityId const entity1, MolecularEntityId const entity2)
{
    return m_topology.remove_edge(entity1, entity2);
}

std::optional<index_t> Topology::first_link(MolecularEntityId const entity, MolecularEntityCategory const category) const
{
    if (!m_topology.contains(entity))
    {
        return std::nullopt;
    }

    for (auto const& [other_category, other_entity] : m_topology.adjacency(entity))
    {
        if (other_category == category)
        {
            return other_entity;
        }
    }

    return std::nullopt;
}

std::vector<index_t> Topology::all_links(MolecularEntityId const entity, MolecularEntityCategory const category) const
{
    std::vector<index_t> links;
    if (!m_topology.contains(entity))
    {
        return links;
    }

    for (auto const& [other_category, other_entity] : m_topology.adjacency(entity))
    {
        if (other_category == category)
        {
            links.push_back(other_entity);
        }
    }

    return links;
}

size_t Topology::count_links(MolecularEntityId const entity, MolecularEntityCategory const category) const
{
    if (!m_topology.contains(entity))
    {
        return 0;
    }

    size_t count = 0;
    for (auto const& [other_category, other_entity] : m_topology.adjacency(entity))
    {
        if (other_category == category)
        {
            ++count;
        }
    }

    return count;
}

bool Topology::contains_link(MolecularEntityId const entity1, MolecularEntityId const entity2) const
{
    return m_topology.contains_edge(entity1, entity2);
}

std::size_t Topology::EntityHash::operator()(MolecularEntityId const& entity) const
{
    // Simple hash that should be good enough for our case.
    // Based on Boost's hash_combine modified by Adam Nevraumont (see https://stackoverflow.com/questions/5889238/why-is-xor-the-default-way-to-combine-hashes/27952689#27952689).
    int const first = static_cast<int>(entity.category);
    index_t const second = entity.index;
    if constexpr (sizeof(size_t) >= 8)
    {
        return first ^ (second + 0x517cc1b727220a95 + (first << 6) + (first >> 2));
    }
    else
    {
        return first ^ (second + 0x9e3779b9 + (first << 6) + (first >> 2));
    }
}

} // namespace mol::internal
