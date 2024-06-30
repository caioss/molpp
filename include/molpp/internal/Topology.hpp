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
    struct MolecularEntityId
    {
        MolecularEntityCategory category;
        index_t index;

        friend bool operator==(MolecularEntityId const& lhs, MolecularEntityId const& rhs);
    };

    bool link_entities(MolecularEntityId const entity1, MolecularEntityId const entity2);
    bool link_categories(MolecularEntityCategory const category1, MolecularEntityCategory const category2);
    bool remove_link(MolecularEntityId const entity1, MolecularEntityId const entity2);
    std::optional<index_t> find_link(MolecularEntityId const entity, MolecularEntityCategory const category) const;
    std::vector<index_t> all_links(MolecularEntityId const entity, MolecularEntityCategory const category) const;
    size_t count_links(MolecularEntityId const entity, MolecularEntityCategory const category) const;
    bool contains_link(MolecularEntityId const entity1, MolecularEntityId const entity2) const;
    std::vector<index_t> convert(std::vector<index_t> const& source_indices, MolecularEntityCategory const source_category, MolecularEntityCategory const target_category);

private:
    using CategoryPair = std::pair<MolecularEntityCategory, MolecularEntityCategory>;

    struct EntityHash
    {
        std::size_t operator()(MolecularEntityId const& entity) const;
    };

    struct CategoryHash
    {
        std::size_t operator()(std::pair<MolecularEntityCategory, MolecularEntityCategory> ends) const;
    };

    SimpleGraph<MolecularEntityId, EntityHash> m_topology;
    SimpleGraph<MolecularEntityCategory> m_hierarchy;
    std::unordered_map<CategoryPair, std::vector<MolecularEntityCategory>, CategoryHash> m_conversion_paths;
};

} // namespace mol::internal

#endif // MOLPP_INTERNAL_TOPOLOGY_HPP
