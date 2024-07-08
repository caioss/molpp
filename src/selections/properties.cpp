#include "selections/properties.hpp"
#include <molpp/internal/MolData.hpp>

namespace mol::internal
{

mol::internal::ResidSelection::ResidSelection(NumberSet&& resids)
: m_resids{resids}
{
}

bool ResidSelection::evaluate_atom(index_t const atom_index, MolData const& data) const
{
    Topology const& topology = data.topology();
    std::optional<index_t> const residue_index = topology.find_link({EntityCategory::Atom, atom_index}, EntityCategory::Residue);
    if (!residue_index)
    {
        return false;
    }

    int const resid = data.residues().id(*residue_index);
    return m_resids.has(resid);
}

} // namespace mol::internal
