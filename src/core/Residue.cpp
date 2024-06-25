#include <molpp/Residue.hpp>
#include <molpp/Atom.hpp>
#include <molpp/internal/MolData.hpp>

namespace mol
{

MolecularEntityCategory Residue::category()
{
    return MolecularEntityCategory::Residue;
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

std::string const& Residue::segid() const
{
    return data().residues().segid(index());
}

void Residue::set_segid(std::string const& segid)
{
    data().residues().segid(index()) = segid;
}

std::string const& Residue::chain() const
{
    return data().residues().chain(index());
}

void Residue::set_chain(std::string const& chain)
{
    data().residues().chain(index()) = chain;
}

void Residue::add_atom(index_t atom_index)
{
    internal::Topology& topology = data().topology();

    std::optional<index_t> const old_residue = topology.find_link({MolecularEntityCategory::Atom, atom_index}, Residue::category());
    if (old_residue)
    {
        topology.remove_link({MolecularEntityCategory::Atom, atom_index}, {Residue::category(), *old_residue});
    }

    topology.link_entities({MolecularEntityCategory::Atom, atom_index}, {Residue::category(), index()});
}

void Residue::add_atom(Atom const& atom)
{
    add_atom(atom.index());
}

size_t Residue::size() const
{
    return data().topology().count_links({Residue::category(), index()}, MolecularEntityCategory::Atom);
}

} // namespace mol
