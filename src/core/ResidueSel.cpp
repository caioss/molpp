#include "molpp/ResidueSel.hpp"
#include <molpp/internal/MolData.hpp>
#include "molpp/internal/SelIndex.hpp"

using namespace mol;
using namespace mol::internal;

size_t ResidueSel::data_size(internal::MolData const& data)
{
    return data.residues().size();
}

std::vector<index_t> ResidueSel::atom_indices() const
{
    ResidueData const &residues = data()->residues();
    size_t num_atoms = 0;
    for (auto const res : indices())
    {
        num_atoms += residues.size(res);
    }

    std::vector<index_t> atoms;
    atoms.reserve(num_atoms);
    for (auto const res : indices())
    {
        for (auto const index : residues.indices(res))
        {
            atoms.push_back(index);
        }
    }

    return atoms;
}

size_t ResidueSel::atom_index(size_t const atom_index, internal::MolData const& data)
{
    return data.atoms().residue(atom_index);
}
