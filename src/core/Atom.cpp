#include <molpp/internal/MolData.hpp>
#include <molpp/Atom.hpp>
#include <molpp/Residue.hpp>
#include <molpp/Error.hpp>

#include <ranges>

namespace mol
{

MolecularEntityCategory Atom::category()
{
    return MolecularEntityCategory::Atom;
}

std::optional<Residue> Atom::residue()
{
    std::optional<index_t> const index = residue_index();
    if (!index)
    {
        return std::nullopt;
    }

    return Residue(*index, frame(), data());
}

std::optional<index_t> Atom::residue_index() const
{
    internal::Topology const& topology = data().topology();
    return topology.find_link({Atom::category(), index()}, Residue::category());
}

int Atom::atomic_number() const
{
    return data().atoms().atomic_number(index());
}

void Atom::set_atomic_number(int const atomic)
{
    data().atoms().atomic_number(index()) = atomic;
}

float Atom::occupancy() const
{
    return data().atoms().occupancy(index());
}

void Atom::set_occupancy(float const occupancy)
{
    data().atoms().occupancy(index()) = occupancy;
}

float Atom::temperature_factor() const
{
    return data().atoms().temperature_factor(index());
}

void Atom::set_temperature_factor(float const tempfactor)
{
    data().atoms().temperature_factor(index()) = tempfactor;
}

float Atom::mass() const
{
    return data().atoms().mass(index());
}

void Atom::set_mass(float const mass)
{
    data().atoms().mass(index()) = mass;
}

float Atom::charge() const
{
    return data().atoms().charge(index());
}

void Atom::set_charge(float const charge)
{
    data().atoms().charge(index()) = charge;
}

float Atom::radius() const
{
    return data().atoms().radius(index());
}

void Atom::set_radius(float const radius)
{
    data().atoms().radius(index()) = radius;
}

std::string const& Atom::name() const
{
    return data().atoms().name(index());
}

void Atom::set_name(std::string const& name)
{
    data().atoms().name(index()) = name;
}

std::string const& Atom::type() const
{
    return data().atoms().type(index());
}

void Atom::set_type(std::string const& type)
{
    data().atoms().type(index()) = type;
}

std::string const& Atom::alternate_location() const
{
    return data().atoms().alternate_location(index());
}

void Atom::set_alternate_location(std::string const& altloc)
{
    data().atoms().alternate_location(index()) = altloc;
}

std::string const& mol::Atom::insertion_code() const
{
    return data().atoms().insertion_code(index());
}

void mol::Atom::set_insertion_code(std::string const& insertion_code)
{
    data().atoms().insertion_code(index()) = insertion_code;
}

std::shared_ptr<Bond> Atom::add_bond(index_t const bonded_to)
{
    if (bonded_to == index())
    {
        throw mol::Error("Atoms can't have bonds to themselves");
    }
    if (bonded_to >= data().size<Atom>())
    {
        throw mol::Error("Out of bounds index: " + std::to_string(bonded_to));
    }
    return data().bonds().add_bond(index(), bonded_to);
}

std::shared_ptr<Bond> Atom::add_bond(Atom const& bonded_to)
{
    return add_bond(bonded_to.index());
}

std::shared_ptr<Bond> Atom::bond(index_t const other)
{
    return data().bonds().bond(index(), other);
}

std::shared_ptr<Bond> Atom::bond(Atom const& other)
{
    return bond(other.index());
}

std::vector<std::shared_ptr<Bond>> mol::Atom::bonds()
{
    std::ranges::single_view indices{index()};
    return data().bonds().bonds(indices.begin(), indices.end());
}

Positions3::ColXpr mol::Atom::position()
{
    if (!frame())
    {
        throw mol::Error("Invalid frame");
    }
    return data().trajectory().timestep(*frame()).coords().col(index());
}

Positions3::ConstColXpr mol::Atom::position() const
{
    if (!frame())
    {
        throw mol::Error("Invalid frame");
    }
    return data().trajectory().timestep(*frame()).coords().col(index());
}

} // namespace mol
