#include "molpp/ResidueSel.hpp"
#include <molpp/internal/MolData.hpp>

using namespace mol;
using namespace mol::internal;

ResidueSel::indices_type ResidueSel::as_atom_indices() const
{
    ResidueData const& residues = data().residues();
    size_t num_atoms = 0;
    for (auto const res : indices())
    {
        num_atoms += residues.size(res);
    }

    ResidueSel::indices_type atoms;
    atoms.reserve(num_atoms);
    for (auto const res : indices())
    {
        for (auto const index : residues.atom_indices(res))
        {
            atoms.push_back(index);
        }
    }

    return atoms;
}
