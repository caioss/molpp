#ifndef ATOMDATA_HPP
#define ATOMDATA_HPP

#include "Timestep.hpp"
#include "BondGraph.hpp"
#include "ResidueData.hpp"
#include "AtomProperties.hpp"
#include <memory>
#include <vector>

namespace mol::internal {

class AtomData : public std::enable_shared_from_this<AtomData>
{
public:
    AtomData() = delete;
    static std::shared_ptr<AtomData> create(size_t const num_atoms);

    size_t size() const { return m_num_atoms; };
    AtomProperties &properties() { return m_properties; }
    BondGraph &bonds() { return m_bonds; }
    ResidueData &residues() { return m_residues; }

    size_t num_frames() { return m_timestep.size(); }
    Timestep &timestep(size_t const index) { return m_timestep[index]; }
    void add_timestep(Timestep &&ts);

private:
    AtomData(size_t const num_atoms);

    size_t m_num_atoms;
    AtomProperties m_properties;
    std::vector<Timestep> m_timestep;
    BondGraph m_bonds;
    ResidueData m_residues;
};

} // namespace mol::internal

#endif // ATOMDATA_HPP
