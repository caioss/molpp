#ifndef MOLDATA_HPP
#define MOLDATA_HPP

#include "core/Timestep.hpp"
#include "core/AtomData.hpp"
#include "core/BondData.hpp"
#include "core/ResidueData.hpp"
#include <memory>
#include <vector>

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

    size_t num_frames() { return m_timestep.size(); }
    Timestep &timestep(size_t const index) { return m_timestep[index]; }
    void add_timestep(Timestep &&ts);

private:
    MolData(size_t const num_atoms);

    size_t m_num_atoms;
    AtomData m_properties;
    std::vector<Timestep> m_timestep;
    BondData m_bonds;
    ResidueData m_residues;
};

} // namespace mol::internal

#endif // MOLDATA_HPP
