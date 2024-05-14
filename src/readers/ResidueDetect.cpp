#include "readers/ResidueDetect.hpp"
#include <molpp/internal/MolData.hpp>

namespace mol::internal
{

ResidueDetect::ResidueDetect()
: m_iterator{m_residues.end()}
{}

index_t ResidueDetect::register_atom(int const resid, std::string const& resname, std::string const& segid, std::string const chain)
{
    if (m_iterator == m_residues.end() || m_current.resid != resid || m_current.resname != resname || m_current.segid != segid || m_current.chain != chain)
    {
        index_t const index = m_residues.size();
        m_current = {index, 0, resid, resname, segid, chain};

        std::tuple const key{resid, resname, segid, chain};
        m_iterator = m_residues.insert(std::pair(key, m_current)).first;
    }

    ResidueDetect::ResidueInfo& residue = m_iterator->second;
    residue.count++;
    return residue.index;
}

void ResidueDetect::update_residue_data(MolData& mol_data) const
{
    ResidueData& residues_data = mol_data.residues();
    residues_data.resize(m_residues.size());
    for (auto const& item : m_residues)
    {
        ResidueDetect::ResidueInfo const& residue = item.second;
        residues_data.clear_and_reserve(residue.index, residue.count);
        residues_data.set(residue.index, residue.resid, residue.resname, residue.segid, residue.chain);
    }

    for (index_t index = 0; index < mol_data.size<Atom>(); ++index)
    {
        index_t const residue_idx = mol_data.atoms().residue(index);
        residues_data.add_atom(residue_idx, index);
    }
}

} // namespace mol::internal
