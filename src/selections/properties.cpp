#include "selections/properties.hpp"
#include "selections/SelectionStack.hpp"
#include "selections/SelectionIndices.hpp"
#include <molpp/Common.hpp>
#include <molpp/internal/MolData.hpp>

namespace mol::internal
{

void PropSelection::evaluate(SelectionStack& stack, MolData const& data, Frame /*frame*/) const
{
    SelectionIndices indices = stack.pop_indices();

    for (index_t atom_index : *(indices.available))
    {
        if (selected(atom_index, data))
        {
            indices.selected->insert(atom_index);
        }
    }
}

bool ResidSelection::selected(index_t atom_index, MolData const& data) const
{
    Topology const& topology = data.topology();
    std::optional<index_t> const residue_index = topology.find_link({EntityCategory::Atom, atom_index}, EntityCategory::Residue);
    if (!residue_index)
    {
        return false;
    }

    int const resid = data.residues().id(*residue_index);
    return has(resid);
}

} // namespace mol::internal
