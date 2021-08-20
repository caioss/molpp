#include <molpp/internal/BaseSel.hpp>
#include <molpp/AtomSel.hpp>
#include "core/MolData.hpp"
#include <numeric>

using namespace mol;
using namespace mol::internal;

BaseSel::BaseSel(size_t const max_size, std::shared_ptr<MolData> data)
: m_data(data),
  m_selected(max_size, true),
  m_indices(max_size)
{
    std::iota(m_indices.begin(), m_indices.end(), 0);
    init_frame();
}

BaseSel::BaseSel(size_t const max_size, std::vector<size_t> const &indices, std::shared_ptr<MolData> data)
: m_data(data),
  m_selected(max_size, false),
  m_indices(indices)
{
    update_indices(max_size);
    init_frame();
}

BaseSel::BaseSel(size_t const max_size, std::vector<size_t> &&indices, std::shared_ptr<MolData> data)
: m_data(data),
  m_selected(max_size, false),
  m_indices(0)
{
    std::swap(m_indices, indices);
    update_indices(max_size);
    init_frame();
}

void BaseSel::set_frame(std::optional<size_t> frame)
{
    if (frame && frame >= m_data->num_frames())
    {
        throw mol::MolError("Out of bounds frame: " + std::to_string(*frame));
    }
    m_frame = frame;
}

BaseSel::coords_type BaseSel::coords(std::vector<size_t> const &atom_indices)
{
    if (!m_frame)
    {
        throw mol::MolError("Invalid frame");
    }
    return m_data->timestep(*m_frame).coords()(Eigen::all, atom_indices);
}

std::vector<size_t> BaseSel::bonded(std::vector<size_t> const &atom_indices) const
{
    return m_data->bonds().bonded(atom_indices.begin(), atom_indices.end());
}

std::vector<std::shared_ptr<mol::Bond>> BaseSel::bonds(std::vector<size_t> const &atom_indices)
{
    return m_data->bonds().bonds(atom_indices.begin(), atom_indices.end());
}

std::shared_ptr<AtomSel> BaseSel::atoms(std::vector<size_t> &&atom_indices)
{
    auto sel = std::make_shared<AtomSel>(AtomSel::from_atom_indices(std::forward<std::vector<size_t>&&>(atom_indices), m_data));
    sel->set_frame(m_frame);
    return sel;
}

void BaseSel::update_indices(size_t const total_size)
{
    std::sort(m_indices.begin(), m_indices.end());
    for (size_t index : m_indices)
    {
        if (index >= total_size)
        {
            throw mol::MolError("Out of bounds selection index: " + std::to_string(index));
        }

        m_selected[index] = true;
    }
}

void BaseSel::init_frame()
{
    if (m_data->num_frames())
    {
        m_frame = 0;
    }
}
