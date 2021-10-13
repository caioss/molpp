#include "tools/SpatialSearch.hpp"
#include <cmath>

using namespace mol::internal;

SpatialSearch::SpatialSearch(Eigen::Ref<point_data_t> const &points, float const cell_size)
: m_cell_size(cell_size),
  m_points{points}
{
    update();
}

std::vector<std::pair<SpatialSearch::index_t, SpatialSearch::index_t>> SpatialSearch::pairs(float const distance) const
{
    index_t const num_layers = floor(distance / m_cell_size) + 1;
    float const distance2 = distance * distance;
    std::vector<std::pair<index_t, index_t>> pairs_list;

    for (index_t cell_z = 1; cell_z < m_grid_size(2) - 1; cell_z++)
    for (index_t cell_y = 1; cell_y < m_grid_size(1) - 1; cell_y++)
    for (index_t cell_x = 1; cell_x < m_grid_size(0) - 1; cell_x++)
    {
    {
    {
        cell_index_t const current_index{cell_x, cell_y, cell_z};
        cell_t const &current_data = m_cells[m_strides * current_index];

        for (index_t dz = -num_layers; dz <= num_layers; dz++)
        for (index_t dy = -num_layers; dy <= num_layers; dy++)
        for (index_t dx = -num_layers; dx <= num_layers; dx++)
        {
        {
        {
            cell_index_t const diff{dx, dy, dz};
            size_t const offset = m_strides * (current_index + diff);
            if (offset >= m_cells.size())
            {
                // big num_layers values may cause overflow
                // of the grid indices
                continue;
            }
            find_cells_pairs(pairs_list, current_data, m_cells[offset], distance2);
        }
        }
        }
    }
    }
    }

    return pairs_list;
}

std::vector<SpatialSearch::index_t> SpatialSearch::query(index_t const query, float const distance) const
{
    float const distance2 = distance * distance;
    index_t const num_layers = floor(distance / m_cell_size) + 1;
    point_t point = m_points.col(query);
    cell_index_t const point_cell = clamped_index(point);
    std::vector<index_t> result;

    for (index_t dz = -num_layers; dz <= num_layers; dz++)
    for (index_t dy = -num_layers; dy <= num_layers; dy++)
    for (index_t dx = -num_layers; dx <= num_layers; dx++)
    {
    {
    {
        cell_index_t const diff{dx, dy, dz};
        size_t const offset = m_strides * (point_cell + diff);
        if (offset >= m_cells.size())
        {
            // big num_layers values may cause overflow
            // of the grid indices
            continue;
        }

        for (index_t const clamped_index : m_cells[offset])
        {
            if ((m_points.col(clamped_index) - point).squaredNorm() <= distance2)
            {
                result.push_back(clamped_index);
            }
        }
    }
    }
    }

    return result;
}

void SpatialSearch::update()
{
    point_t min_coords = m_points.rowwise().minCoeff();
    point_t max_coords = m_points.rowwise().maxCoeff();

    // Avoid huge grids
    float const max_edge = max_coords.maxCoeff() - min_coords.minCoeff();
    if (max_edge / m_cell_size > 100)
    {
        m_cell_size = max_edge / 100;
    }

    // Extra layers are added to help queries on the boundaries
    m_origin = min_coords.array() - m_cell_size;
    m_max_clamp = index(max_coords);
    m_grid_size = m_max_clamp.array() + 2;
    m_strides = {1, m_grid_size(0), m_grid_size(0) * m_grid_size(1)};
    m_cells.resize(m_grid_size.prod());

    // Pre-reserving saves considerable time at the
    // expense of more memory use. Don't reserve on
    // the boundary cells.
    float const density = m_points.cols() / m_max_clamp.prod();
    for (index_t cell_z = 1; cell_z < m_grid_size(2) - 1; cell_z++)
    for (index_t cell_y = 1; cell_y < m_grid_size(1) - 1; cell_y++)
    for (index_t cell_x = 1; cell_x < m_grid_size(0) - 1; cell_x++)
    {
    {
    {
        size_t const offset = m_strides * cell_index_t{cell_x, cell_y, cell_z};
        m_cells[offset].reserve(density);
    }
    }
    }

    // Categorize all points
    for (index_t i = 0; i < m_points.cols(); i++)
    {
        size_t const offset = m_strides * clamped_index(m_points.col(i));
        m_cells[offset].push_back(i);
    }
}

SpatialSearch::cell_index_t SpatialSearch::index(Eigen::Ref<point_t> const &point) const
{
    return ((point - m_origin) / m_cell_size).array().floor().cast<index_t>();
}

SpatialSearch::cell_index_t SpatialSearch::clamped_index(Eigen::Ref<point_t> const &point) const
{
    return index(std::forward<Eigen::Ref<point_t> const>(point)).cwiseMax(cell_index_t{1, 1, 1}).cwiseMin(m_max_clamp);
}

void SpatialSearch::find_cells_pairs(std::vector<std::pair<index_t, index_t>> &pairs_list, cell_t const &current, cell_t const &neighbor, float const distance2) const
{
    for (index_t const i : current)
    for (index_t const j : neighbor)
    {
    {
        if (i <= j)
        {
            // Avoid duplicates
            continue;
        }

        if ((m_points.col(i) - m_points.col(j)).squaredNorm() <= distance2)
        {
            pairs_list.push_back(std::pair(i, j));
        }
    }
    }
}
