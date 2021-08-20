#include "core/MolData.hpp"
#include <molpp/ResidueSel.hpp>

using namespace mol;

ResidueSel ResidueSel::from_atom_indices(std::vector<size_t> &&indices, std::shared_ptr<mol::internal::MolData> data)
{
    std::unordered_set<size_t> residues;
    for (auto const &index : indices)
    {
        residues.insert(data->properties().residue(index));
    }

    return ResidueSel(std::vector<size_t>(residues.begin(), residues.end()), data);
}

#if 0
std::vector<size_t> ResidueSel::from_atom_indices(std::vector<size_t> &&indices)
{
    std::unordered_set<size_t> residues;
    for (auto const &index : indices)
    {
        residues.insert(data->properties().residue(index));
    }

    return {residues.begin(), residues.end()};
}
#endif

std::vector<size_t> ResidueSel::atom_indices()
{
    size_t num_atoms = 0;
    for (auto &res : indices())
    {
        num_atoms += data()->residues().size(res);
    }

    std::vector<size_t> atoms;
    atoms.reserve(num_atoms);
    for (auto const &res : indices())
    {
        for (auto &index : data()->residues().indices(res))
        {
            atoms.push_back(index);
        }
    }

    return atoms;
}

size_t ResidueSel::max_size(std::shared_ptr<mol::internal::MolData> data)
{
    return data->residues().size();
}
