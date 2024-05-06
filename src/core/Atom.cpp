#include <molpp/internal/MolData.hpp>
#include <molpp/Atom.hpp>
#include <molpp/Residue.hpp>
#include <molpp/MolError.hpp>

#include <ranges>

using namespace mol;

int Atom::resid() const
{
    return data()->residues().residue_id(residue_id());
}

Residue Atom::residue()
{
    return Residue(residue_id(), frame(), data());
}

index_t Atom::residue_id() const
{
    return data()->atoms().residue(index());
}

int Atom::atomic_number() const
{
    return data()->atoms().atomic_number(index());
}

void Atom::set_atomic_number(int const atomic)
{
    data()->atoms().atomic_number(index()) = atomic;
}

float Atom::occupancy() const
{
    return data()->atoms().occupancy(index());
}

void Atom::set_occupancy(float const occupancy)
{
    data()->atoms().occupancy(index()) = occupancy;
}

float Atom::temperature_factor() const
{
    return data()->atoms().temperature_factor(index());
}

void Atom::set_temperature_factor(float const tempfactor)
{
    data()->atoms().temperature_factor(index()) = tempfactor;
}

float Atom::mass() const
{
    return data()->atoms().mass(index());
}

void Atom::set_mass(float const mass)
{
    data()->atoms().mass(index()) = mass;
}

float Atom::charge() const
{
    return data()->atoms().charge(index());
}

void Atom::set_charge(float const charge)
{
    data()->atoms().charge(index()) = charge;
}

float Atom::radius() const
{
    return data()->atoms().radius(index());
}

void Atom::set_radius(float const radius)
{
    data()->atoms().radius(index()) = radius;
}

std::string const& Atom::name() const
{
    return data()->atoms().name(index());
}

void Atom::set_name(std::string const& name)
{
    data()->atoms().name(index()) = name;
}

std::string const& Atom::type() const
{
    return data()->atoms().type(index());
}

void Atom::set_type(std::string const& type)
{
    data()->atoms().type(index()) = type;
}

std::string const& Atom::residue_name() const
{
    return data()->residues().residue_name(residue_id());
}

std::string const& Atom::segid() const
{
    return data()->residues().segid(residue_id());
}

std::string const& Atom::chain() const
{
    return data()->residues().chain(residue_id());
}

std::string const& Atom::alternate_location() const
{
    return data()->atoms().alternate_location(index());
}

void Atom::set_alternate_location(std::string const& altloc)
{
    data()->atoms().alternate_location(index()) = altloc;
}

std::string const& mol::Atom::insertion_code() const
{
    return data()->atoms().insertion_code(index());
}

void mol::Atom::set_insertion_code(std::string const& insertion_code)
{
    data()->atoms().insertion_code(index()) = insertion_code;
}

std::shared_ptr<Bond> Atom::add_bond(index_t const bonded_to)
{
    if (bonded_to == index())
    {
        throw mol::MolError("Atoms can't have bonds to themselves");
    }
    if (bonded_to >= data()->size<Atom>())
    {
        throw mol::MolError("Out of bounds index: " + std::to_string(bonded_to));
    }
    return data()->bonds().add_bond(index(), bonded_to);
}

std::shared_ptr<Bond> Atom::add_bond(Atom const& bonded_to)
{
    return add_bond(bonded_to.index());
}

std::shared_ptr<Bond> Atom::bond(index_t const other)
{
    return data()->bonds().bond(index(), other);
}

std::shared_ptr<Bond> Atom::bond(Atom const& other)
{
    return bond(other.index());
}

std::vector<index_t> Atom::as_atom_indices() const
{
    return {index()};
}

std::vector<std::shared_ptr<Bond>> mol::Atom::bonds()
{
    std::ranges::single_view indices{index()};
    return data()->bonds().bonds(indices.begin(), indices.end());
}
