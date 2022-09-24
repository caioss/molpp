#include <molpp/internal/BaseAtomAggregate.hpp>
#include <molpp/MolError.hpp>
#include "core/MolData.hpp"

using namespace mol;
using namespace mol::internal;

BaseAtomAggregate::coords_type BaseAtomAggregate::coords(std::vector<index_t> &&atom_indices)
{
    if (!m_frame)
    {
        throw mol::MolError("Invalid frame");
    }
    return m_data->trajectory().timestep(*m_frame).coords()(Eigen::all, std::forward<std::vector<index_t>>(atom_indices));
}

std::vector<std::shared_ptr<Bond>> BaseAtomAggregate::bonds(std::vector<index_t> const &atom_indices)
{
    return m_data->bonds().bonds(atom_indices.begin(), atom_indices.end());
}
