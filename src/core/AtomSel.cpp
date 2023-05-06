#include <molpp/AtomSel.hpp>
#include "molpp/internal/SelIndex.hpp"
#include "core/MolData.hpp"

using namespace mol;

size_t AtomSel::data_size(internal::MolData const& data)
{
    return data.size();
}

std::vector<index_t> const& AtomSel::atom_indices() const
{
    return indices();
}
