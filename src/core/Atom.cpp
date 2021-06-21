#include "Atom.hpp"
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

int Atom::resid() const
{
    return m_data->m_resid[m_index];
}

void Atom::set_resid(int const &resid)
{
    m_data->m_resid[m_index] = resid;
}

int Atom::atomic() const
{
    return m_data->m_atomic[m_index];
}

void Atom::set_atomic(int const &atomic)
{
    m_data->m_atomic[m_index] = atomic;
}

float Atom::occupancy() const
{
    return m_data->m_occupancy[m_index];
}

void Atom::set_occupancy(float const &occupancy)
{
    m_data->m_occupancy[m_index] = occupancy;
}

float Atom::tempfactor() const
{
    return m_data->m_tempfactor[m_index];
}

void Atom::set_tempfactor(float const &tempfactor)
{
    m_data->m_tempfactor[m_index] = tempfactor;
}

float Atom::mass() const
{
    return m_data->m_mass[m_index];
}

void Atom::set_mass(float const &mass)
{
    m_data->m_mass[m_index] = mass;
}

float Atom::charge() const
{
    return m_data->m_charge[m_index];
}

void Atom::set_charge(float const &charge)
{
    m_data->m_charge[m_index] = charge;
}

float Atom::radius() const
{
    return m_data->m_radius[m_index];
}

void Atom::set_radius(float const &radius)
{
    m_data->m_radius[m_index] = radius;
}

std::string Atom::name() const
{
    return m_data->m_name[m_index];
}

void Atom::set_name(std::string const &name) {
    m_data->m_name[m_index] = name;
}

std::string Atom::type() const
{
    return m_data->m_type[m_index];
}

void Atom::set_type(std::string const &type)
{
    m_data->m_type[m_index] = type;
}

std::string Atom::resname() const
{
    return m_data->m_resname[m_index];
}

void Atom::set_resname(std::string const &resname)
{
    m_data->m_resname[m_index] = resname;
}

std::string Atom::segid() const
{
    return m_data->m_segid[m_index];
}

void Atom::set_segid(std::string const &segid)
{
    m_data->m_segid[m_index] = segid;
}

std::string Atom::chain() const
{
    return m_data->m_chain[m_index];
}

void Atom::set_chain(std::string const &chain)
{
    m_data->m_chain[m_index] = chain;
}

std::string Atom::altloc() const
{
    return m_data->m_altloc[m_index];
}

void Atom::set_altloc(std::string const &altloc)
{
    m_data->m_altloc[m_index] = altloc;
}

Eigen::Ref<Pos3> Atom::coords()
{
    return m_data->m_timestep[m_frame].coords().col(m_index);
}
