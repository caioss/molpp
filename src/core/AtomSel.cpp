#include "core/MolData.hpp"
#include <molpp/AtomSel.hpp>

using namespace mol;

AtomSel::coords_type AtomSel::coords()
{
    return m_data->timestep(frame()).coords()(Eigen::all, indices());
}

std::shared_ptr<AtomSel> AtomSel::bonded()
{
    return std::make_shared<AtomSel>(m_data->size(), m_data->bonds().bonded(indices().begin(), indices().end()), m_data);
}

std::vector<std::shared_ptr<Bond>> AtomSel::bonds()
{
    return m_data->bonds().bonds(indices().begin(), indices().end());
}
