#include "readers/ResidueDetect.hpp"
#include <molpp/MolppCore.hpp>
#include <molpp/internal/MolData.hpp>

namespace mol::internal
{

ResidueDetect::ResidueDetect(MolData& data)
: m_data{data}
{
}

void ResidueDetect::register_atom(index_t const atom_index, int const resid, std::string const& resname, std::string const& segid, std::string const& chain)
{
    auto chain_iter = m_chains.find(chain);
    if (chain_iter == m_chains.end())
    {
        chain_iter = m_chains.emplace(chain, m_chains.size()).first;
        m_chain_name.push_back(chain);
    }

    auto segment_iter = m_segments.find(segid);
    if (segment_iter == m_segments.end())
    {
        segment_iter = m_segments.emplace(segid, m_segments.size()).first;
        m_segment_name.push_back(segid);
    }

    index_t const chain_index = chain_iter->second;
    index_t const segment_index = segment_iter->second;

    ResidueKey query{resid, chain_index, segment_index, resname};
    auto residue_iter = m_residues.find(query);
    if (residue_iter == m_residues.end())
    {
        residue_iter = m_residues.emplace(query, m_residues.size()).first;
    }

    index_t const residue_index = residue_iter->second;
    Topology& topology = m_data.topology();
    topology.link_entities({MolecularEntityCategory::Residue, residue_index}, {MolecularEntityCategory::Atom, atom_index});
}

void ResidueDetect::update_residue_data(MolData& mol_data) const
{
    // Update the hierarchy
    Topology& topology = mol_data.topology();
    if (!m_residues.empty())
    {
        topology.link_categories(MolecularEntityCategory::Residue, MolecularEntityCategory::Atom);
    }

    ResidueData& residues_data = mol_data.residues();
    residues_data.resize(m_residues.size());
    for (auto const& item : m_residues)
    {
        ResidueKey const& residue = item.first;
        residues_data.set(item.second, residue.resid, residue.resname, m_segment_name[residue.segment], m_chain_name[residue.chain]);
    }
}

bool operator<(ResidueDetect::ResidueKey const& lhs, ResidueDetect::ResidueKey const& rhs)
{
    return std::tie(lhs.resid, lhs.resname, lhs.segment, lhs.chain) < std::tie(rhs.resid, rhs.resname, rhs.segment, rhs.chain);
}

} // namespace mol::internal
