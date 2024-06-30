#include <molpp/internal/Topology.hpp>
#include "tools/algorithms.hpp"

size_t pair_hash(size_t const first, size_t const second)
{
    // Simple hash that should be good enough for our case.
    // Based on Boost's hash_combine modified by Adam Nevraumont (see https://stackoverflow.com/questions/5889238/why-is-xor-the-default-way-to-combine-hashes/27952689#27952689).
    if constexpr (sizeof(size_t) >= 8)
    {
        return first ^ (second + 0x517cc1b727220a95 + (first << 6) + (first >> 2));
    }
    else
    {
        return first ^ (second + 0x9e3779b9 + (first << 6) + (first >> 2));
    }
}

namespace mol::internal
{

bool operator==(Topology::EntityId const& lhs, Topology::EntityId const& rhs)
{
    return lhs.category == rhs.category && lhs.index == rhs.index;
}

bool Topology::link_entities(EntityId const entity1, EntityId const entity2)
{
    return m_topology.add_edge(entity1, entity2, true);
}

bool Topology::link_categories(EntityCategory const category1, EntityCategory const category2)
{
    return m_hierarchy.add_edge(category1, category2, true);
}

bool Topology::remove_link(EntityId const entity1, EntityId const entity2)
{
    return m_topology.remove_edge(entity1, entity2);
}

std::optional<index_t> Topology::find_link(EntityId const entity, EntityCategory const category) const
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

std::vector<index_t> Topology::all_links(EntityId const entity, EntityCategory const category) const
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

size_t Topology::count_links(EntityId const entity, EntityCategory const category) const
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

bool Topology::contains_link(EntityId const entity1, EntityId const entity2) const
{
    return m_topology.contains_edge(entity1, entity2);
}

std::vector<index_t> Topology::convert(std::vector<index_t> const& source_indices, EntityCategory const source_category, EntityCategory const target_category)
{
    if (source_category == target_category)
    {
        return source_indices;
    }

    std::vector<EntityCategory>& path = m_conversion_paths[{source_category, target_category}];

    if (path.empty())
    {
        BreadthFirstTraversal bfs(m_hierarchy);
        bool const result = bfs.run(source_category, [=](EntityCategory const category) {
            return category == target_category;
        }, [=](auto) {
            return true;
        });

        if (!result)
        {
            return {};
        }

        // Build both paths
        std::vector<EntityCategory>& inverse_path = m_conversion_paths[{target_category, source_category}];
        auto& parents = bfs.parent_map();

        EntityCategory current = parents.at(target_category);
        while (current != source_category)
        {
            inverse_path.push_back(current);
            current = parents.at(current);
        }

        path.resize(inverse_path.size() + 1);
        std::copy(inverse_path.rbegin(), inverse_path.rend(), path.begin());

        // Add last step
        path.back() = target_category;
        inverse_path.push_back(source_category);
    }

    EntityCategory current_category = source_category;
    std::vector<index_t> indices = source_indices;
    std::vector<index_t> new_indices;
    for (EntityCategory const next_category : path)
    {
        new_indices.clear();
        for (index_t const index : indices)
        {
            std::vector<index_t> const links = all_links({current_category, index}, next_category);
            new_indices.insert(new_indices.end(), links.begin(), links.end());
        }

        std::swap(indices, new_indices);
        current_category = next_category;
    }

    return indices;
}

std::size_t Topology::EntityHash::operator()(EntityId const& entity) const
{
    return pair_hash(static_cast<size_t>(entity.category), entity.index);
}

std::size_t Topology::CategoryHash::operator()(std::pair<EntityCategory, EntityCategory> ends) const
{
    return pair_hash(static_cast<size_t>(ends.first), static_cast<size_t>(ends.second));
}

} // namespace mol::internal
