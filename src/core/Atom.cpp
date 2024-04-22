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

int Atom::atomic() const
{
    return data()->atoms().atomic(index());
}

void Atom::set_atomic(int const &atomic)
{
    data()->atoms().atomic(index()) = atomic;
}

float Atom::occupancy() const
{
    return data()->atoms().occupancy(index());
}

void Atom::set_occupancy(float const &occupancy)
{
    data()->atoms().occupancy(index()) = occupancy;
}

float Atom::tempfactor() const
{
    return data()->atoms().tempfactor(index());
}

void Atom::set_tempfactor(float const &tempfactor)
{
    data()->atoms().tempfactor(index()) = tempfactor;
}

float Atom::mass() const
{
    return data()->atoms().mass(index());
}

void Atom::set_mass(float const &mass)
{
    data()->atoms().mass(index()) = mass;
}

float Atom::charge() const
{
    return data()->atoms().charge(index());
}

void Atom::set_charge(float const &charge)
{
    data()->atoms().charge(index()) = charge;
}

float Atom::radius() const
{
    return data()->atoms().radius(index());
}

void Atom::set_radius(float const &radius)
{
    data()->atoms().radius(index()) = radius;
}

std::string Atom::name() const
{
    return data()->atoms().name(index());
}

void Atom::set_name(std::string const &name) {
    data()->atoms().name(index()) = name;
}

std::string Atom::type() const
{
    return data()->atoms().type(index());
}

void Atom::set_type(std::string const &type)
{
    data()->atoms().type(index()) = type;
}

std::string Atom::resname() const
{
    return data()->residues().resname(residue_id());
}

std::string Atom::segid() const
{
    return data()->residues().segid(residue_id());
}

std::string Atom::chain() const
{
    return data()->residues().chain(residue_id());
}

std::string Atom::altloc() const
{
    return data()->atoms().altloc(index());
}

void Atom::set_altloc(std::string const &altloc)
{
    data()->atoms().altloc(index()) = altloc;
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
    return index() < data()->atoms().size();
}
