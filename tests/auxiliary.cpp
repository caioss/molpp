#include "auxiliary.hpp"

using namespace mol;
using namespace mol::internal;

std::shared_ptr<MolData> aux_moldata()
{
    size_t const num_atoms { 3 };
    auto data = MolData::create(num_atoms);
    data->residues().resize(3);
    for (size_t i = 0; i < num_atoms; i++)
    {
        data->properties().residue(i) = i;
        data->residues().add_atom(i, i);
    }

    data->bonds().add_bond(0, 2);

    data->add_timestep(Timestep(num_atoms));
    data->timestep(0).coords() << 1, 2, 3, 4, 5, 6, 7, 8, 9;
    return data;
}
