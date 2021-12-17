#include "core/MolData.hpp"
#include <molpp/ResidueSel.hpp>

using namespace mol;
using namespace mol::internal;

std::vector<size_t> ResidueSel::atom_indices() const
{
    ResidueData const &residues = cdata()->residues();
    size_t num_atoms = 0;
    for (auto const res : indices())
    {
        num_atoms += residues.size(res);
    }

    std::vector<size_t> atoms;
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

std::vector<size_t> ResidueSel::from_atom_indices(std::vector<size_t> &&atom_indices, std::shared_ptr<mol::internal::MolData> const data)
{
    AtomData &properties = data->properties();
    std::unordered_set<size_t> residues;
    for (auto const index : atom_indices)
    {
        residues.insert(properties.residue(index));
    }

    return {residues.begin(), residues.end()};
}

size_t ResidueSel::max_size(std::shared_ptr<mol::internal::MolData> const data)
{
    return data->residues().size();
}
