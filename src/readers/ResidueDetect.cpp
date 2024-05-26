#include "readers/ResidueDetect.hpp"
#include <molpp/internal/MolData.hpp>

namespace mol::internal
{

index_t ResidueDetect::register_atom(int const resid, std::string const& resname, std::string const& segid, std::string const& chain)
{
    auto chain_iter = m_chains.find(chain);
    if (chain_iter == m_chains.end())
    {
        chain_iter = m_chains.emplace(chain, EntityInfo{m_chains.size(), 0}).first;
        m_chain_name.push_back(chain);
    }

    auto segment_iter = m_segments.find(segid);
    if (segment_iter == m_segments.end())
    {
        segment_iter = m_segments.emplace(segid, EntityInfo{m_segments.size(), 0}).first;
        m_segment_name.push_back(segid);
    }

    EntityInfo& chain_info = chain_iter->second;
    EntityInfo& segment_info = segment_iter->second;

    ResidueKey query{resid, chain_info.index, segment_info.index, resname};
    auto residue_iter = m_residues.find(query);
    if (residue_iter == m_residues.end())
    {
        residue_iter = m_residues.emplace(query, EntityInfo{m_residues.size(), 0}).first;
        chain_info.size++;
        segment_info.size++;
    }

    residue_iter->second.size++;
    return residue_iter->second.index;
}

void ResidueDetect::update_residue_data(MolData& mol_data) const
{
    ResidueData& residues_data = mol_data.residues();
    residues_data.resize(m_residues.size());
    for (auto const& item : m_residues)
    {
        ResidueKey const& residue = item.first;
        EntityInfo const& info = item.second;
        residues_data.clear_and_reserve(info.index, info.size);
        residues_data.set(info.index, residue.resid, residue.resname, m_segment_name[residue.segment], m_chain_name[residue.chain]);
    }

    mol::internal::AtomData& atom_data = mol_data.atoms();
    for (index_t atom_index = 0; atom_index < mol_data.size<Atom>(); ++atom_index)
    {
        std::optional<index_t> const residue_index = atom_data.residue(atom_index);
        if (!residue_index)
        {
            continue;
        }
        residues_data.add_atom(*residue_index, atom_index);
    }
}

bool operator<(ResidueDetect::ResidueKey const& lhs, ResidueDetect::ResidueKey const& rhs)
{
    return std::tie(lhs.resid, lhs.resname, lhs.segment, lhs.chain) < std::tie(rhs.resid, rhs.resname, rhs.segment, rhs.chain);
}

} // namespace mol::internal
