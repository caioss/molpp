#include <molpp/internal/BaseSel.hpp>
#include <molpp/AtomSel.hpp>
#include <molpp/ResidueSel.hpp>
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
  m_selected(max_size, false)
{
    // Fill m_selected
    size_t count = 0;
    for (size_t index : indices)
    {
        if (index >= max_size)
        {
            throw mol::MolError("Out of bounds selection index: " + std::to_string(index));
        }

        if (!m_selected[index])
        {
            m_selected[index] = true;
            count++;
        }
    }

    // Fill m_indices with sorted and non-duplicated values
    m_indices.reserve(count);
    for (size_t i = 0; i < max_size; i++)
    {
        if (m_selected[i])
        {
            m_indices.push_back(i);
        }
    }

    init_frame();
}

void BaseSel::set_frame(std::optional<size_t> frame)
{
    if (frame && frame >= m_data->trajectory().num_frames())
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
    return m_data->trajectory().timestep(*m_frame).coords()(Eigen::all, atom_indices);
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
    auto sel = std::make_shared<AtomSel>(AtomSel::from_atom_indices(std::forward<std::vector<size_t>>(atom_indices), m_data));
    sel->set_frame(m_frame);
    return sel;
}

std::shared_ptr<ResidueSel> BaseSel::residues(std::vector<size_t> &&atom_indices)
{
    auto sel = std::make_shared<ResidueSel>(ResidueSel::from_atom_indices(std::forward<std::vector<size_t>>(atom_indices), m_data));
    sel->set_frame(m_frame);
    return sel;
}

void BaseSel::init_frame()
{
    if (m_data->trajectory().num_frames())
    {
        m_frame = 0;
    }
}
