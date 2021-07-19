#include "Atom.hpp"
#include "AtomSel.hpp"
#include "MolError.hpp"
#include "core/AtomData.hpp"

using namespace mol;

Atom::Atom(size_t const index, size_t const frame, std::shared_ptr<internal::AtomData> data)
: m_index { index },
  m_frame { frame },
  m_data { data }
{
}

bool Atom::operator==(Atom const &other) const
{
    return m_data == other.m_data && m_index == other.m_index && m_frame == other.m_frame;
}

size_t Atom::index() const
{
    return m_index;
}

size_t Atom::frame() const
{
    return m_frame;
}

int Atom::resid() const
{
    return m_data->properties().resid(m_index);
}

void Atom::set_resid(int const &resid)
{
    m_data->properties().resid(m_index) = resid;
}

int Atom::atomic() const
{
    return m_data->properties().atomic(m_index);
}

void Atom::set_atomic(int const &atomic)
{
    m_data->properties().atomic(m_index) = atomic;
}

float Atom::occupancy() const
{
    return m_data->properties().occupancy(m_index);
}

void Atom::set_occupancy(float const &occupancy)
{
    m_data->properties().occupancy(m_index) = occupancy;
}

float Atom::tempfactor() const
{
    return m_data->properties().tempfactor(m_index);
}

void Atom::set_tempfactor(float const &tempfactor)
{
    m_data->properties().tempfactor(m_index) = tempfactor;
}

float Atom::mass() const
{
    return m_data->properties().mass(m_index);
}

void Atom::set_mass(float const &mass)
{
    m_data->properties().mass(m_index) = mass;
}

float Atom::charge() const
{
    return m_data->properties().charge(m_index);
}

void Atom::set_charge(float const &charge)
{
    m_data->properties().charge(m_index) = charge;
}

float Atom::radius() const
{
    return m_data->properties().radius(m_index);
}

void Atom::set_radius(float const &radius)
{
    m_data->properties().radius(m_index) = radius;
}

std::string Atom::name() const
{
    return m_data->properties().name(m_index);
}

void Atom::set_name(std::string const &name) {
    m_data->properties().name(m_index) = name;
}

std::string Atom::type() const
{
    return m_data->properties().type(m_index);
}

void Atom::set_type(std::string const &type)
{
    m_data->properties().type(m_index) = type;
}

std::string Atom::resname() const
{
    return m_data->properties().resname(m_index);
}

void Atom::set_resname(std::string const &resname)
{
    m_data->properties().resname(m_index) = resname;
}

std::string Atom::segid() const
{
    return m_data->properties().segid(m_index);
}

void Atom::set_segid(std::string const &segid)
{
    m_data->properties().segid(m_index) = segid;
}

std::string Atom::chain() const
{
    return m_data->properties().chain(m_index);
}

void Atom::set_chain(std::string const &chain)
{
    m_data->properties().chain(m_index) = chain;
}

std::string Atom::altloc() const
{
    return m_data->properties().altloc(m_index);
}

void Atom::set_altloc(std::string const &altloc)
{
    m_data->properties().altloc(m_index) = altloc;
}

std::shared_ptr<Bond> Atom::add_bond(size_t const bonded_to)
{
    if (bonded_to == m_index)
    {
        throw mol::MolError("Atoms can't have bonds to themselves");
    }
    if (bonded_to >= m_data->size())
    {
        throw mol::MolError("Out of bounds index: " + std::to_string(bonded_to));
    }
    return m_data->bonds().add_bond(m_index, bonded_to);
}

std::shared_ptr<Bond> Atom::bond(size_t const other)
{
    return m_data->bonds().bond(m_index, other);
}

std::vector<std::shared_ptr<Bond>> Atom::bonds()
{
    return m_data->bonds().bonds(m_index);
}

std::shared_ptr<AtomSel> Atom::bonded()
{
    return std::make_shared<AtomSel>(m_data->size(), m_data->bonds().bonded(m_index), m_data);
}

Eigen::Ref<Pos3> Atom::coords()
{
    return m_data->timestep(m_frame).coords().col(m_index);
}
