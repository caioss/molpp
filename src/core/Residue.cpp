#include "core/MolData.hpp"
#include <molpp/Residue.hpp>
#include <molpp/AtomSel.hpp>

using namespace mol;

int Residue::resid() const
{
    return cdata()->residues().resid(index());
}

void Residue::set_resid(int const &resid)
{
    data()->residues().resid(index()) = resid;
}

std::string Residue::resname() const
{
    return cdata()->residues().resname(index());
}

void Residue::set_resname(std::string const &resname)
{
    data()->residues().resname(index()) = resname;
}

std::string Residue::segid() const
{
    return cdata()->residues().segid(index());
}

void Residue::set_segid(std::string const &segid)
{
    data()->residues().segid(index()) = segid;
}

std::string Residue::chain() const
{
    return cdata()->residues().chain(index());
}

void Residue::set_chain(std::string const &chain)
{
    data()->residues().chain(index()) = chain;
}

void Residue::add_atom(size_t atom_index)
{
    mol::internal::AtomData &properties = data()->properties();
    mol::internal::ResidueData &residues = data()->residues();
    size_t const old_res = properties.residue(atom_index);
    residues.remove_atom(old_res, atom_index);
    residues.add_atom(index(), atom_index);
    properties.residue(atom_index) = index();
}

void Residue::add_atom(Atom const &atom)
{
    add_atom(atom.index());
}

size_t Residue::size() const
{
    return cdata()->residues().size(index());
}

std::vector<size_t> Residue::atom_indices() const
{
    auto const indices = cdata()->residues().indices(index());
    return {indices.begin(), indices.end()};
}
