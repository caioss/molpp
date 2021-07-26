#ifndef ATOMPROPERTIES_HPP
#define ATOMPROPERTIES_HPP

#include <vector>
#include <string>

namespace mol::internal {

class AtomProperties
{
public:
    AtomProperties() = delete;

    AtomProperties(size_t const num_atoms)
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
    size_t &residue(size_t const index) { return m_residue[index]; }
    size_t &atomic(size_t const index) { return m_atomic[index]; }
    float &occupancy(size_t const index) { return m_occupancy[index]; }
    float &tempfactor(size_t const index) { return m_tempfactor[index]; }
    float &mass(size_t const index) { return m_mass[index]; }
    float &charge(size_t const index) { return m_charge[index]; }
    float &radius(size_t const index) { return m_radius[index]; }
    std::string &name(size_t const index) { return m_name[index]; }
    std::string &type(size_t const index) { return m_type[index]; }
    std::string &altloc(size_t const index) { return m_altloc[index]; }

private:
    size_t m_num_atoms;
    std::vector<size_t> m_residue;
    std::vector<size_t> m_atomic;
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

#endif // ATOMPROPERTIES_HPP
