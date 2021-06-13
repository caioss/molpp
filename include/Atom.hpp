#ifndef ATOM_HPP
#define ATOM_HPP

#include "core/AtomData.hpp"
#include <Eigen/Dense>
#include <memory>

namespace mol {

class Atom
{
public:
    Atom(size_t const index, size_t const frame,
         std::shared_ptr<internal::AtomData> data)
    : m_index { index },
      m_frame { frame },
      m_data { data }
    {}

    bool operator==(Atom const &other) const { return m_data == other.m_data && m_index == other.m_index && m_frame == other.m_frame; }

    size_t index() { return m_index; }

    int resid() const { return m_data->m_resid[m_index]; }
    void set_resid(int const &resid) { m_data->m_resid[m_index] = resid; }

    int atomic() const { return m_data->m_atomic[m_index]; }
    void set_atomic(int const &atomic) { m_data->m_atomic[m_index] = atomic; }

    float occupancy() const { return m_data->m_occupancy[m_index]; }
    void set_occupancy(float const &occupancy) { m_data->m_occupancy[m_index] = occupancy; }

    float tempfactor() const { return m_data->m_tempfactor[m_index]; }
    void set_tempfactor(float const &tempfactor) { m_data->m_tempfactor[m_index] = tempfactor; }

    float mass() const { return m_data->m_mass[m_index]; }
    void set_mass(float const &mass) { m_data->m_mass[m_index] = mass; }

    float charge() const { return m_data->m_charge[m_index]; }
    void set_charge(float const &charge) { m_data->m_charge[m_index] = charge; }

    float radius() const { return m_data->m_radius[m_index]; }
    void set_radius(float const &radius) { m_data->m_radius[m_index] = radius; }

    std::string name() const { return m_data->m_name[m_index]; }
    void set_name(std::string const &name) { m_data->m_name[m_index] = name; }

    std::string type() const { return m_data->m_type[m_index]; }
    void set_type(std::string const &type) { m_data->m_type[m_index] = type; }

    std::string resname() const { return m_data->m_resname[m_index]; }
    void set_resname(std::string const &resname) { m_data->m_resname[m_index] = resname; }

    std::string segid() const { return m_data->m_segid[m_index]; }
    void set_segid(std::string const &segid) { m_data->m_segid[m_index] = segid; }

    std::string chain() const { return m_data->m_chain[m_index]; }
    void set_chain(std::string const &chain) { m_data->m_chain[m_index] = chain; }

    std::string altloc() const { return m_data->m_altloc[m_index]; }
    void set_altloc(std::string const &altloc) { m_data->m_altloc[m_index] = altloc; }

    Eigen::Ref<Eigen::Vector3f> coords() { return m_data->m_timestep[m_frame].coords().col(m_index); }

private:
    size_t m_index;
    size_t m_frame;
    std::shared_ptr<internal::AtomData> m_data;
};

} // namespace mol

#endif // ATOM_HPP
