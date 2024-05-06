#include <molpp/internal/MolData.hpp>
#include <molpp/Residue.hpp>
#include <molpp/AtomSel.hpp>

using namespace mol;

int Residue::residue_id() const
{
    return data()->residues().residue_id(index());
}

void Residue::set_residue_id(int const resid)
{
    data()->residues().residue_id(index()) = resid;
}

std::string const& Residue::residue_name() const
{
    return data()->residues().residue_name(index());
}

void Residue::set_residue_name(std::string const& resname)
{
    data()->residues().residue_name(index()) = resname;
}

std::string const& Residue::segid() const
{
    return data()->residues().segid(index());
}

void Residue::set_segid(std::string const& segid)
{
    data()->residues().segid(index()) = segid;
}

std::string const& Residue::chain() const
{
    return data()->residues().chain(index());
}

void Residue::set_chain(std::string const& chain)
{
    data()->residues().chain(index()) = chain;
}

void Residue::add_atom(index_t atom_index)
{
    mol::internal::AtomData& properties = data()->atoms();
    mol::internal::ResidueData& residues = data()->residues();
    index_t const old_res = properties.residue(atom_index);
    residues.remove_atom(old_res, atom_index);
    residues.add_atom(index(), atom_index);
    properties.residue(atom_index) = index();
}

void Residue::add_atom(Atom const& atom)
{
    add_atom(atom.index());
}

size_t Residue::size() const
{
    return data()->residues().size(index());
}

std::vector<index_t> Residue::as_atom_indices() const
{
    auto const indices = data()->residues().indices(index());
    return {indices.begin(), indices.end()};
}
