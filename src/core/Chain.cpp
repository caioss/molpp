#include <molpp/Chain.hpp>
#include <molpp/Residue.hpp>
#include <molpp/internal/MolData.hpp>

namespace mol
{

MolecularEntityCategory Chain::category()
{
    return MolecularEntityCategory::Chain;
}

size_t Chain::size() const
{
    return data().topology().count_links({category(), index()}, MolecularEntityCategory::Residue);
}

std::string const& Chain::name() const
{
    return data().chains().name(index());
}

void Chain::set_name(std::string const& resname)
{
    data().chains().name(index()) = resname;
}

void Chain::add_residue(index_t residue_index)
{
    internal::Topology& topology = data().topology();
    internal::Topology::MolecularEntityId const residue_id{MolecularEntityCategory::Residue, residue_index};

    std::optional<index_t> const old_chain = topology.find_link(residue_id, category());
    if (old_chain)
    {
        topology.remove_link(residue_id, {category(), *old_chain});
    }

    topology.link_entities(residue_id, {category(), index()});
}

void Chain::add_residue(Residue const& residue)
{
    add_residue(residue.index());
}

} // namespace mol
