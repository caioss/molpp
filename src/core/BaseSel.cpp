#include <molpp/internal/BaseSel.hpp>
#include <molpp/AtomSel.hpp>
#include <molpp/ResidueSel.hpp>
#include "core/MolData.hpp"
#include <numeric>

using namespace mol;
using namespace mol::internal;

BaseSel::BaseSel(SelIndex&& indices, std::shared_ptr<MolData> data)
: m_data(data),
  m_index{indices}
{
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

BaseSel::coords_type BaseSel::coords(std::vector<size_t> &&atom_indices)
{
    if (!m_frame)
    {
        throw mol::MolError("Invalid frame");
    }
    return m_data->trajectory().timestep(*m_frame).coords()(Eigen::all, std::forward<std::vector<size_t>>(atom_indices));
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
    auto sel = std::make_shared<AtomSel>(AtomSel::from_atom_indices(std::forward<std::vector<size_t>>(atom_indices), m_data), m_data);
    sel->set_frame(m_frame);
    return sel;
}

std::shared_ptr<ResidueSel> BaseSel::residues(std::vector<size_t> &&atom_indices)
{
    auto sel = std::make_shared<ResidueSel>(ResidueSel::from_atom_indices(std::forward<std::vector<size_t>>(atom_indices), m_data), m_data);
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
