#include "selections/properties.hpp"
#include "selections/SelectionStack.hpp"
#include <molpp/Common.hpp>
#include <molpp/internal/MolData.hpp>

namespace mol::internal
{

void PropSelection::evaluate(SelectionStack& stack, MolData const& data, Frame /*frame*/) const
{
    SelectionFlags flags = stack.pop_flags();

    for (index_t atom_index : *(flags.mask))
    {
        if (selected(atom_index, data))
        {
            flags.selected->insert(atom_index);
        }
    }
}

bool ResidSelection::selected(index_t atom_index, MolData const& data) const
{
    Topology const& topology = data.topology();
    std::optional<index_t> const residue_index = topology.find_link({MolecularEntityCategory::Atom, atom_index}, MolecularEntityCategory::Residue);
    if (!residue_index)
    {
        return false;
    }

    int const resid = data.residues().id(*residue_index);
    return has(resid);
}

} // namespace mol::internal
