#include <molpp/internal/MolData.hpp>
#include <molpp/Atom.hpp>
#include <molpp/Residue.hpp>
#include <molpp/MolError.hpp>

using namespace mol;

int Atom::resid() const
{
    return data()->residues().resid(residue_id());
}

Residue Atom::residue()
{
    return Residue(residue_id(), frame(), data());
}

index_t Atom::residue_id() const
{
    return data()->atoms().residue(index());
}

std::shared_ptr<Bond> Atom::add_bond(index_t const bonded_to)
{
    if (bonded_to == index())
    {
        throw mol::MolError("Atoms can't have bonds to themselves");
    }
    if (bonded_to >= data()->size())
    {
        throw mol::MolError("Out of bounds index: " + std::to_string(bonded_to));
    }
    return data()->bonds().add_bond(index(), bonded_to);
}

std::shared_ptr<Bond> Atom::add_bond(Atom const &bonded_to)
{
    return add_bond(bonded_to.index());
}

std::shared_ptr<Bond> Atom::bond(index_t const other)
{
    return data()->bonds().bond(index(), other);
}

std::shared_ptr<Bond> Atom::bond(Atom const &other)
{
    return bond(other.index());
}

std::vector<index_t> Atom::atom_indices() const
{
    return {index()};
}

bool Atom::validate_index() const
{
    return index() < data()->properties().size<Atom>();
}
