#ifndef MOLDATA_HPP
#define MOLDATA_HPP

#include "core/AtomData.hpp"
#include "core/BondData.hpp"
#include "core/ResidueData.hpp"
#include "core/TrajData.hpp"
#include <memory>

namespace mol::internal {

class MolData : public std::enable_shared_from_this<MolData>
{
public:
    MolData() = delete;
    static std::shared_ptr<MolData> create(size_t const num_atoms);

    size_t size() const { return m_num_atoms; };
    AtomData &properties() { return m_properties; }
    BondData &bonds() { return m_bonds; }
    ResidueData &residues() { return m_residues; }
    TrajData &trajectory() { return m_trajectory; }

private:
    MolData(size_t const num_atoms);

    size_t m_num_atoms;
    AtomData m_properties;
    BondData m_bonds;
    ResidueData m_residues;
    TrajData m_trajectory;
};

} // namespace mol::internal

#endif // MOLDATA_HPP
