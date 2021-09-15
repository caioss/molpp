#include <molpp/internal/BaseAtomAggregate.hpp>
#include <molpp/AtomSel.hpp>
#include "core/MolData.hpp"

using namespace mol;
using namespace mol::internal;

BaseAtomAggregate::BaseAtomAggregate(size_t const index, std::optional<size_t> const frame,std::shared_ptr<internal::MolData> data)
: m_index { index },
  m_frame(frame),
  m_data(data)
{}

bool BaseAtomAggregate::operator==(BaseAtomAggregate const &other) const
{
    return m_data == other.m_data
           && m_index == other.m_index
           && m_frame == other.m_frame;
}

BaseAtomAggregate::coords_type BaseAtomAggregate::coords(std::vector<size_t> &&atom_indices)
{
    if (!m_frame)
    {
        throw mol::MolError("Invalid frame");
    }
    return m_data->trajectory().timestep(*m_frame).coords()(Eigen::all, std::forward<std::vector<size_t>>(atom_indices));
}

std::vector<std::shared_ptr<Bond>> BaseAtomAggregate::bonds(std::vector<size_t> const &atom_indices)
{
    return m_data->bonds().bonds(atom_indices.begin(), atom_indices.end());
}

std::shared_ptr<AtomSel> BaseAtomAggregate::bonded(std::vector<size_t> const &atom_indices)
{
    auto &&bonded = m_data->bonds().bonded(atom_indices.begin(), atom_indices.end());
    auto sel = std::make_shared<AtomSel>(std::forward<std::vector<size_t>>(bonded), m_data);
    sel->set_frame(frame());
    return sel;
}

std::shared_ptr<AtomSel> BaseAtomAggregate::atoms(std::vector<size_t> &&atom_indices)
{
    auto sel = std::make_shared<AtomSel>(std::forward<std::vector<size_t>>(atom_indices), m_data);
    sel->set_frame(frame());
    return sel;
}
