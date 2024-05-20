#include <molpp/AtomSel.hpp>
#include <molpp/internal/MolData.hpp>

namespace mol
{

AtomSel::position_type AtomSel::positions()
{
    return data().trajectory().timestep(*frame()).coords()(Eigen::all, indices());
}

AtomSel::const_position_type AtomSel::positions() const
{
    return data().trajectory().timestep(*frame()).coords()(Eigen::all, indices());
}

AtomSel::indices_type const& AtomSel::as_atom_indices() const
{
    return indices();
}

AtomSel AtomSel::bonded()
{
    std::vector<index_t> bonded_atoms = data().bonds().bonded(indices().begin(), indices().end());
    AtomSel sel(bonded_atoms, data());
    sel.set_frame(frame());
    return sel;
}

std::vector<std::shared_ptr<Bond>> AtomSel::bonds()
{
    return data().bonds().bonds(indices().begin(), indices().end());
}

} // namespace mol
