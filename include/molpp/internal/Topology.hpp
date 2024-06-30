#ifndef MOLPP_INTERNAL_TOPOLOGY_HPP
#define MOLPP_INTERNAL_TOPOLOGY_HPP

#include <molpp/Common.hpp>
#include <molpp/tools/SimpleGraph.hpp>

#include <vector>
#include <utility>
#include <optional>
#include <unordered_map>

namespace mol::internal
{

class Topology
{
public:
    struct EntityId
    {
        EntityCategory category;
        index_t index;

        friend bool operator==(EntityId const& lhs, EntityId const& rhs);
    };

    bool link_entities(EntityId const entity1, EntityId const entity2);
    bool link_categories(EntityCategory const category1, EntityCategory const category2);
    bool remove_link(EntityId const entity1, EntityId const entity2);
    std::optional<index_t> find_link(EntityId const entity, EntityCategory const category) const;
    std::vector<index_t> all_links(EntityId const entity, EntityCategory const category) const;
    size_t count_links(EntityId const entity, EntityCategory const category) const;
    bool contains_link(EntityId const entity1, EntityId const entity2) const;
    std::vector<index_t> convert(std::vector<index_t> const& source_indices, EntityCategory const source_category, EntityCategory const target_category);

private:
    using CategoryPair = std::pair<EntityCategory, EntityCategory>;

    struct EntityHash
    {
        std::size_t operator()(EntityId const& entity) const;
    };

    struct CategoryHash
    {
        std::size_t operator()(std::pair<EntityCategory, EntityCategory> ends) const;
    };

    SimpleGraph<EntityId, EntityHash> m_topology;
    SimpleGraph<EntityCategory> m_hierarchy;
    std::unordered_map<CategoryPair, std::vector<EntityCategory>, CategoryHash> m_conversion_paths;
};

} // namespace mol::internal

#endif // MOLPP_INTERNAL_TOPOLOGY_HPP
