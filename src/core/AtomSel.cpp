#include <molpp/AtomSel.hpp>
#include <molpp/internal/SelIndex.hpp>
#include <molpp/internal/MolData.hpp>

using namespace mol;

size_t AtomSel::data_size(internal::MolData const& data)
{
    return data.size();
}

std::vector<index_t> const& AtomSel::atom_indices() const
{
    return indices();
}

mol::AtomSel mol::AtomSel::bonded()
{
    std::vector<index_t> bonded_atoms = data()->bonds().bonded(indices().begin(), indices().end());
    AtomSel sel(from_atom_indices(bonded_atoms, *data()), data());
    sel.set_frame(frame());
    return sel;
}

std::vector<std::shared_ptr<mol::Bond>> mol::AtomSel::bonds()
{
    return data()->bonds().bonds(indices().begin(), indices().end());
}
