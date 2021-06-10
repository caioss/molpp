#ifndef ATOMDATA_HPP
#define ATOMDATA_HPP

#include "Timestep.hpp"
#include <memory>
#include <vector>
#include <string>

namespace mol {

class Atom;

namespace internal {

class AtomData : public std::enable_shared_from_this<AtomData>
{
public:
    AtomData() = delete;
    static std::shared_ptr<AtomData> create(size_t const num_atoms);

    size_t size() const { return m_num_atoms; };
    mol::Atom index(size_t const index, size_t const frame);

    size_t num_frames() { return m_timestep.size(); }
    Timestep &timestep(size_t const index) { return m_timestep[index]; }
    void add_timestep(Timestep &&ts);

private:
    AtomData(size_t const num_atoms);

    size_t m_num_atoms;
    std::vector<int> m_resid;
    std::vector<int> m_residue;
    std::vector<int> m_atomic;
    std::vector<float> m_occupancy;
    std::vector<float> m_tempfactor;
    std::vector<float> m_mass;
    std::vector<float> m_charge;
    std::vector<float> m_radius;
    std::vector<std::string> m_name;
    std::vector<std::string> m_type;
    std::vector<std::string> m_resname;
    std::vector<std::string> m_segid;
    std::vector<std::string> m_chain;
    std::vector<std::string> m_altloc;
    std::vector<Timestep> m_timestep;

    friend class mol::Atom;
};

} // namespace internal
} // namespace mol

#endif // ATOMDATA_HPP
