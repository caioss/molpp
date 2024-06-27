#include "readers/ResidueDetect.hpp"
#include <molpp/MolppCore.hpp>
#include <molpp/internal/MolData.hpp>

#include <limits>

namespace mol::internal
{

// Use the maximum value of index_t to avoid overlapping with valid indices
constexpr index_t const INVALID_INDEX = std::numeric_limits<index_t>::max();

ResidueDetect::ResidueDetect(MolData& data)
: m_data{data}
{
    // Register the empty chain and segment to allow computing further indices
    m_chains[""] = INVALID_INDEX;
    m_segments[""] = INVALID_INDEX;
}

void ResidueDetect::register_atom(index_t const atom_index, int const resid, std::string const& resname, std::string const& segid, std::string const& chain)
{
    auto segment_iter = m_segments.find(segid);
    if (segment_iter == m_segments.end())
    {
        segment_iter = m_segments.emplace(segid, m_segments.size() - 1).first;
        m_segment_name.push_back(segid);
    }

    index_t const chain_index = register_chain(chain);
    index_t const segment_index = segment_iter->second;

    Topology& topology = m_data.topology();
    ResidueKey query{resid, chain_index, segment_index, resname};
    auto residue_iter = m_residues.find(query);
    if (residue_iter == m_residues.end())
    {
        residue_iter = m_residues.emplace(query, m_residues.size()).first;

        if (!chain.empty())
        {
            topology.link_entities({MolecularEntityCategory::Chain, chain_index}, {MolecularEntityCategory::Residue, residue_iter->second});
        }
    }

    index_t const residue_index = residue_iter->second;
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

    if (m_chains.size() > 1)
    {
        topology.link_categories(MolecularEntityCategory::Chain, MolecularEntityCategory::Residue);
    }

    ChainData& chain_data = mol_data.chains();
    chain_data.resize(m_chains.size() - 1);
    for (auto const& item : m_chains)
    {
        if (item.first.empty())
        {
            continue;
        }

        chain_data.name(item.second) = item.first;
    }

    ResidueData& residues_data = mol_data.residues();
    residues_data.resize(m_residues.size());
    for (auto const& item : m_residues)
    {
        ResidueKey const& residue = item.first;
        index_t const index = item.second;
        residues_data.id(index) = residue.resid;
        residues_data.name(index) = residue.resname;
    }
}

index_t ResidueDetect::register_chain(std::string const& chain)
{
    auto chain_iter = m_chains.find(chain);
    if (chain_iter == m_chains.end())
    {
        index_t const index = m_chains.size() - 1; // Remove the empty one
        chain_iter = m_chains.emplace(chain, index).first;
        return index;
    }

    return chain_iter->second;
}

bool operator<(ResidueDetect::ResidueKey const& lhs, ResidueDetect::ResidueKey const& rhs)
{
    return std::tie(lhs.resid, lhs.resname, lhs.segment, lhs.chain) < std::tie(rhs.resid, rhs.resname, rhs.segment, rhs.chain);
}

} // namespace mol::internal
