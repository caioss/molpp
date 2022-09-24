#ifndef ATOMDATA_HPP
#define ATOMDATA_HPP

#include <molpp/MolppCore.hpp>
#include <vector>
#include <string>

namespace mol::internal {

class AtomData
{
public:
    AtomData() = delete;

    AtomData(size_t const num_atoms)
    : m_num_atoms { num_atoms },
      m_residue(num_atoms, -1),
      m_atomic(num_atoms, 0),
      m_occupancy(num_atoms, 0),
      m_tempfactor(num_atoms, 0),
      m_mass(num_atoms, 0),
      m_charge(num_atoms, 0),
      m_radius(num_atoms, 0),
      m_name(num_atoms),
      m_type(num_atoms),
      m_altloc(num_atoms)
    {}

    size_t size() const { return m_num_atoms; }
    index_t &residue(index_t const index) { return m_residue[index]; }
    index_t const &residue(index_t const index) const { return m_residue[index]; }
    index_t &atomic(index_t const index) { return m_atomic[index]; }
    index_t const &atomic(index_t const index) const { return m_atomic[index]; }
    float &occupancy(index_t const index) { return m_occupancy[index]; }
    float const &occupancy(index_t const index) const { return m_occupancy[index]; }
    float &tempfactor(index_t const index) { return m_tempfactor[index]; }
    float const &tempfactor(index_t const index) const { return m_tempfactor[index]; }
    float &mass(index_t const index) { return m_mass[index]; }
    float const &mass(index_t const index) const { return m_mass[index]; }
    float &charge(index_t const index) { return m_charge[index]; }
    float const &charge(index_t const index) const { return m_charge[index]; }
    float &radius(index_t const index) { return m_radius[index]; }
    float const &radius(index_t const index) const { return m_radius[index]; }
    std::string &name(index_t const index) { return m_name[index]; }
    std::string const &name(index_t const index) const { return m_name[index]; }
    std::string &type(index_t const index) { return m_type[index]; }
    std::string const &type(index_t const index) const { return m_type[index]; }
    std::string &altloc(index_t const index) { return m_altloc[index]; }
    std::string const &altloc(index_t const index) const { return m_altloc[index]; }

private:
    size_t m_num_atoms;
    std::vector<index_t> m_residue;
    std::vector<index_t> m_atomic;
    std::vector<float> m_occupancy;
    std::vector<float> m_tempfactor;
    std::vector<float> m_mass;
    std::vector<float> m_charge;
    std::vector<float> m_radius;
    std::vector<std::string> m_name;
    std::vector<std::string> m_type;
    std::vector<std::string> m_altloc;
};

} // namespace mol::internal

#endif // ATOMDATA_HPP
