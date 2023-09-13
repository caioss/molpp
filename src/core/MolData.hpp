#ifndef MOLDATA_HPP
#define MOLDATA_HPP

#include <molpp/Trajectory.hpp>
#include "core/AtomData.hpp"
#include "core/BondData.hpp"
#include "core/ResidueData.hpp"
#include "core/PropertyContainer.hpp"

namespace mol::internal {

class MolData
{
public:
    MolData(size_t const num_atoms);

    size_t size() const { return m_num_atoms; };
    AtomData &atoms() { return m_properties_old; }
    AtomData const &atoms() const { return m_properties_old; }
    BondData &bonds() { return m_bonds; }
    BondData const& bonds() const { return m_bonds; }
    ResidueData &residues() { return m_residues; }
    ResidueData const &residues() const { return m_residues; }
    Trajectory &trajectory() { return m_trajectory; }
    Trajectory const& trajectory() const { return m_trajectory; }

    PropertyContainer& properties()
    {
        return m_properties;
    }

    PropertyContainer const& properties() const
    {
        return m_properties;
    }

private:
    size_t m_num_atoms;
    PropertyContainer m_properties;
    AtomData m_properties_old;
    BondData m_bonds;
    ResidueData m_residues;
    Trajectory m_trajectory;
};

} // namespace mol::internal

#endif // MOLDATA_HPP
