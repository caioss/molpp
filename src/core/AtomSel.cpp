#include <molpp/AtomSel.hpp>
#include "core/MolData.hpp"

using namespace mol;

std::vector<size_t> AtomSel::atom_indices() const
{
    return indices();
}

std::vector<size_t> AtomSel::from_atom_indices(std::vector<size_t> &&atom_indices, std::shared_ptr<mol::internal::MolData> const /*data*/)
{
    return std::forward<std::vector<size_t>>(atom_indices);
}

size_t AtomSel::max_size(std::shared_ptr<mol::internal::MolData> const data)
{
    return data->size();
}
