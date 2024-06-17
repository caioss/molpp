#include "molpp/ResidueSel.hpp"
#include <molpp/internal/MolData.hpp>

namespace mol
{

ResidueSel::indices_type ResidueSel::as_atom_indices() const
{
    internal::ResidueData const& residues = data().residues();
    internal::Topology const& topology = data().topology();
    ResidueSel::indices_type atoms;
    atoms.reserve(residues.size());
    for (auto const residue_index : indices())
    {
        for (auto const atom_index : topology.all_links({MolecularEntityCategory::Residue, residue_index}, MolecularEntityCategory::Atom))
        {
            atoms.push_back(atom_index);
        }
    }

    return atoms;
}

} // namespace mol
