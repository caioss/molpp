#include <molpp/Residue.hpp>
#include <molpp/Atom.hpp>
#include <molpp/Chain.hpp>
#include <molpp/Segment.hpp>
#include <molpp/internal/MolData.hpp>

namespace mol
{

EntityCategory Residue::category()
{
    return EntityCategory::Residue;
}

std::optional<Chain> Residue::chain()
{
    std::optional<index_t> const index = chain_index();
    if (!index)
    {
        return std::nullopt;
    }

    return Chain(*index, frame(), data());
}

std::optional<index_t> Residue::chain_index() const
{
    internal::Topology const& topology = data().topology();
    return topology.find_link({category(), index()}, EntityCategory::Chain);
}

std::optional<Segment> Residue::segment()
{
    std::optional<index_t> const index = segment_index();
    if (!index)
    {
        return std::nullopt;
    }

    return Segment(*index, frame(), data());
}

std::optional<index_t> Residue::segment_index() const
{
    internal::Topology const& topology = data().topology();
    return topology.find_link({category(), index()}, EntityCategory::Segment);
}

int Residue::id() const
{
    return data().residues().id(index());
}

void Residue::set_id(int const resid)
{
    data().residues().id(index()) = resid;
}

std::string const& Residue::name() const
{
    return data().residues().name(index());
}

void Residue::set_name(std::string const& resname)
{
    data().residues().name(index()) = resname;
}

void Residue::add_atom(index_t atom_index)
{
    internal::Topology& topology = data().topology();
    internal::Topology::EntityId const atom_id{EntityCategory::Atom, atom_index};

    std::optional<index_t> const old_residue = topology.find_link(atom_id, category());
    if (old_residue)
    {
        topology.remove_link(atom_id, {category(), *old_residue});
    }

    topology.link_entities(atom_id, {category(), index()});
}

void Residue::add_atom(Atom const& atom)
{
    add_atom(atom.index());
}

size_t Residue::size() const
{
    return data().topology().count_links({category(), index()}, EntityCategory::Atom);
}

} // namespace mol
