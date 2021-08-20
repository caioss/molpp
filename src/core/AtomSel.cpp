#include <molpp/AtomSel.hpp>
#include "core/MolData.hpp"

using namespace mol;

AtomSel AtomSel::from_atom_indices(std::vector<size_t> &&atom_indices, std::shared_ptr<mol::internal::MolData> data)
{
    return AtomSel(std::forward<std::vector<size_t>&&>(atom_indices), data);
}

std::vector<size_t> AtomSel::atom_indices() const
{
    return indices();
}

size_t AtomSel::max_size(std::shared_ptr<mol::internal::MolData> data)
{
    return data->size();
}
