#ifndef MOLPP_INTERNAL_MOLDATA_HPP
#define MOLPP_INTERNAL_MOLDATA_HPP

#include <molpp/internal/Topology.hpp>
#include <molpp/internal/AtomData.hpp>
#include <molpp/internal/BondData.hpp>
#include <molpp/internal/ResidueData.hpp>
#include <molpp/internal/ChainData.hpp>
#include <molpp/Trajectory.hpp>

namespace mol
{
class Atom;
class Residue;
} // namespace mol

namespace mol::internal
{

class MolData
{
public:
    MolData() = delete;
    MolData(size_t const num_atoms);

    template<class Entity>
    size_t size() const;

    Topology& topology()
    {
        return m_topology;
    }

    Topology const& topology() const
    {
        return m_topology;
    }

    AtomData& atoms()
    {
        return m_atoms;
    }

    AtomData const& atoms() const
    {
        return m_atoms;
    }

    BondData& bonds()
    {
        return m_bonds;
    }

    BondData const& bonds() const
    {
        return m_bonds;
    }

    ResidueData& residues()
    {
        return m_residues;
    }

    ResidueData const& residues() const
    {
        return m_residues;
    }

    ChainData& chains()
    {
        return m_chains;
    }

    ChainData const& chains() const
    {
        return m_chains;
    }

    Trajectory& trajectory()
    {
        return m_trajectory;
    }

    Trajectory const& trajectory() const
    {
        return m_trajectory;
    }

private:
    Topology m_topology;
    AtomData m_atoms;
    BondData m_bonds;
    ResidueData m_residues;
    ChainData m_chains;
    Trajectory m_trajectory;
};

template<>
inline size_t MolData::size<mol::Atom>() const
{
    return m_atoms.size();
}

template<>
inline size_t MolData::size<mol::Residue>() const
{
    return m_residues.size();
}

} // namespace mol::internal

#endif // MOLPP_INTERNAL_MOLDATA_HPP
