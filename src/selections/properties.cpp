#include "selections/properties.hpp"
#include "core/MolData.hpp"

using namespace mol;
using namespace mol::internal;

bool ResidSelection::selected(index_t atom_idx, MolData const& data) const
{
    index_t const res_idx = data.properties().residue(atom_idx);
    int const resid = data.residues().resid(res_idx);
    return has(resid);
}
