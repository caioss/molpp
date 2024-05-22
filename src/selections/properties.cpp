#include "selections/properties.hpp"
#include "selections/SelectionStack.hpp"
#include <molpp/internal/MolData.hpp>

using namespace mol;
using namespace mol::internal;

void PropSelection::evaluate(SelectionStack& stack, MolData const& data, Frame /*frame*/) const
{
    SelectionFlags flags = stack.pop_flags();

    for (index_t atom_idx : *(flags.mask))
    {
        if (selected(atom_idx, data))
        {
            flags.selected->insert(atom_idx);
        }
    }
}

bool ResidSelection::selected(index_t atom_idx, MolData const& data) const
{
    std::optional<index_t> const residue_index = data.atoms().residue(atom_idx);
    if (!residue_index)
    {
        return false;
    }

    int const resid = data.residues().residue_id(*residue_index);
    return has(resid);
}
