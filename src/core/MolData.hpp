#ifndef MOLDATA_HPP
#define MOLDATA_HPP

#include "core/AtomData.hpp"
#include "core/BondData.hpp"
#include "core/ResidueData.hpp"
#include <molpp/Trajectory.hpp>

namespace mol::internal {

class MolData
{
public:
    MolData() = delete;
    MolData(size_t const num_atoms);

    size_t size() const { return m_num_atoms; };
    AtomData &atoms() { return m_properties; }
    AtomData const &atoms() const { return m_properties; }
    BondData &bonds() { return m_bonds; }
    BondData const& bonds() const { return m_bonds; }
    ResidueData &residues() { return m_residues; }
    ResidueData const &residues() const { return m_residues; }
    Trajectory &trajectory() { return m_trajectory; }
    Trajectory const& trajectory() const { return m_trajectory; }

private:
    size_t m_num_atoms;
    AtomData m_properties;
    BondData m_bonds;
    ResidueData m_residues;
    Trajectory m_trajectory;
};

} // namespace mol::internal

#endif // MOLDATA_HPP
