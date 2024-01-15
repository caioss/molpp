#include <molpp/AtomSel.hpp>
#include <molpp/internal/SelIndex.hpp>
#include <molpp/internal/MolData.hpp>

using namespace mol;

size_t AtomSel::data_size(internal::MolData const& data)
{
    return data.properties().size<Atom>();
}

std::vector<index_t> const& AtomSel::atom_indices() const
{
    return indices();
}
