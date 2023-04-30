#include <molpp/AtomSel.hpp>
#include "core/MolData.hpp"

using namespace mol;

std::vector<index_t> AtomSel::atom_indices() const
{
    return indices();
}

std::vector<index_t> AtomSel::from_atom_indices(std::vector<index_t>&& atom_indices, internal::MolData const& /*data*/)
{
    return std::forward<std::vector<index_t>>(atom_indices);
}

size_t AtomSel::max_size(internal::MolData const& data)
{
    return data.size();
}
