#ifndef MOLPP_INTERNAL_TOPOLOGY_HPP
#define MOLPP_INTERNAL_TOPOLOGY_HPP

#include <molpp/MolppCore.hpp>
#include <molpp/tools/SimpleGraph.hpp>

#include <optional>

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

    bool add_link(MolecularEntityId const entity1, MolecularEntityId const entity2);
    bool remove_link(MolecularEntityId const entity1, MolecularEntityId const entity2);
    std::optional<index_t> first_link(MolecularEntityId const entity, MolecularEntityCategory const category) const;
    std::vector<index_t> all_links(MolecularEntityId const entity, MolecularEntityCategory const category) const;
    size_t count_links(MolecularEntityId const entity, MolecularEntityCategory const category) const;
    bool contains_link(MolecularEntityId const entity1, MolecularEntityId const entity2) const;

private:
    struct EntityHash
    {
        std::size_t operator()(MolecularEntityId const& entity) const;
    };

    SimpleGraph<MolecularEntityId, EntityHash> m_topology;
    SimpleGraph<MolecularEntityCategory> m_hierarchy;
};

} // namespace mol::internal

#endif // MOLPP_INTERNAL_TOPOLOGY_HPP
