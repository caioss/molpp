#ifndef MOLDATA_HPP
#define MOLDATA_HPP

#include <molpp/internal/AtomData.hpp>
#include <molpp/internal/BondData.hpp>
#include <molpp/internal/ResidueData.hpp>
#include <molpp/internal/PropertyContainer.hpp>

namespace mol::internal {

class MolData
{
public:
    MolData(size_t const num_atoms);

    AtomData &atoms() { return m_properties_old; }
    AtomData const &atoms() const { return m_properties_old; }
    BondData &bonds() { return m_bonds; }
    BondData const& bonds() const { return m_bonds; }
    ResidueData &residues() { return m_residues; }
    ResidueData const &residues() const { return m_residues; }

    PropertyContainer& properties()
    {
        return m_properties;
    }

    PropertyContainer const& properties() const
    {
        return m_properties;
    }

private:
    PropertyContainer m_properties;
    AtomData m_properties_old;
    BondData m_bonds;
    ResidueData m_residues;
};

} // namespace mol::internal

#endif // MOLDATA_HPP
