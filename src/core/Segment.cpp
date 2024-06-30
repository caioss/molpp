#include <molpp/Segment.hpp>
#include <molpp/Residue.hpp>
#include <molpp/internal/MolData.hpp>

namespace mol
{

EntityCategory Segment::category()
{
    return EntityCategory::Segment;
}

size_t Segment::size() const
{
    return data().topology().count_links({category(), index()}, EntityCategory::Residue);
}

std::string const& Segment::name() const
{
    return data().segments().name(index());
}

void Segment::set_name(std::string const& resname)
{
    data().segments().name(index()) = resname;
}

void Segment::add_residue(index_t residue_index)
{
    internal::Topology& topology = data().topology();
    internal::Topology::EntityId const residue_id{EntityCategory::Residue, residue_index};

    std::optional<index_t> const old_segment = topology.find_link(residue_id, category());
    if (old_segment)
    {
        topology.remove_link(residue_id, {category(), *old_segment});
    }

    topology.link_entities(residue_id, {category(), index()});
}

void Segment::add_residue(Residue const& residue)
{
    add_residue(residue.index());
}

} // namespace mol
