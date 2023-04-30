#include <molpp/internal/BaseSel.hpp>
#include <molpp/MolError.hpp>
#include "core/MolData.hpp"
#include <numeric>

using namespace mol;
using namespace mol::internal;

BaseSel::BaseSel(SelIndex&& indices, MolData* data)
: m_data{data},
  m_index{indices}
{
    init_frame();
}

void BaseSel::set_frame(Frame frame)
{
    if (frame && frame >= m_data->trajectory().num_frames())
    {
        throw mol::MolError("Out of bounds frame: " + std::to_string(*frame));
    }
    m_frame = frame;
}

Timestep& BaseSel::timestep()
{
    if (!m_frame)
    {
        throw mol::MolError("Invalid frame");
    }
    return m_data->trajectory().timestep(m_frame.value());
}

BaseSel::coords_type BaseSel::coords(std::vector<index_t> &&atom_indices)
{
    if (!m_frame)
    {
        throw mol::MolError("Invalid frame");
    }
    return m_data->trajectory().timestep(*m_frame).coords()(Eigen::all, std::forward<std::vector<index_t>>(atom_indices));
}

std::vector<index_t> BaseSel::bonded(std::vector<index_t> const &atom_indices) const
{
    return m_data->bonds().bonded(atom_indices.begin(), atom_indices.end());
}

std::vector<std::shared_ptr<mol::Bond>> BaseSel::bonds(std::vector<index_t> const &atom_indices)
{
    return m_data->bonds().bonds(atom_indices.begin(), atom_indices.end());
}

void BaseSel::init_frame()
{
    if (m_data->trajectory().num_frames())
    {
        m_frame = 0;
    }
}
