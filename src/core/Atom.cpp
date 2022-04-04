#include "core/MolData.hpp"
#include <molpp/Atom.hpp>
#include <molpp/Residue.hpp>
#include <molpp/AtomSel.hpp>
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

size_t Atom::residue_id() const
{
    return data()->properties().residue(index());
}

int Atom::atomic() const
{
    return data()->properties().atomic(index());
}

void Atom::set_atomic(int const &atomic)
{
    data()->properties().atomic(index()) = atomic;
}

float Atom::occupancy() const
{
    return data()->properties().occupancy(index());
}

void Atom::set_occupancy(float const &occupancy)
{
    data()->properties().occupancy(index()) = occupancy;
}

float Atom::tempfactor() const
{
    return data()->properties().tempfactor(index());
}

void Atom::set_tempfactor(float const &tempfactor)
{
    data()->properties().tempfactor(index()) = tempfactor;
}

float Atom::mass() const
{
    return data()->properties().mass(index());
}

void Atom::set_mass(float const &mass)
{
    data()->properties().mass(index()) = mass;
}

float Atom::charge() const
{
    return data()->properties().charge(index());
}

void Atom::set_charge(float const &charge)
{
    data()->properties().charge(index()) = charge;
}

float Atom::radius() const
{
    return data()->properties().radius(index());
}

void Atom::set_radius(float const &radius)
{
    data()->properties().radius(index()) = radius;
}

std::string Atom::name() const
{
    return data()->properties().name(index());
}

void Atom::set_name(std::string const &name) {
    data()->properties().name(index()) = name;
}

std::string Atom::type() const
{
    return data()->properties().type(index());
}

void Atom::set_type(std::string const &type)
{
    data()->properties().type(index()) = type;
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
    return data()->properties().altloc(index());
}

void Atom::set_altloc(std::string const &altloc)
{
    data()->properties().altloc(index()) = altloc;
}

std::shared_ptr<Bond> Atom::add_bond(size_t const bonded_to)
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

std::shared_ptr<Bond> Atom::bond(size_t const other)
{
    return data()->bonds().bond(index(), other);
}

std::shared_ptr<Bond> Atom::bond(Atom const &other)
{
    return bond(other.index());
}
