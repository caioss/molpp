#include "selections/properties.hpp"
#include "selections/SelectionStack.hpp"
#include "core/MolData.hpp"

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
    index_t const res_idx = data.atoms().residue(atom_idx);
    int const resid = data.residues().resid(res_idx);
    return has(resid);
}
