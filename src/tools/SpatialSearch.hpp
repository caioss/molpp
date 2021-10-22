#ifndef SPATIALSEARCH_HPP
#define SPATIALSEARCH_HPP

#include <molpp/MolppCore.hpp>
#include <vector>
#include <utility>
#include <cstdlib>

namespace mol::internal {

template <typename T>
class SpatialSearch {
public:
    using index_t = ptrdiff_t;

private:
    using stride_t = Eigen::RowVector3<index_t>;
    using cell_index_t = Eigen::Vector3<index_t>;
    using cell_t = std::vector<index_t>;

public:
    SpatialSearch() = delete;
    SpatialSearch(SpatialSearch const &other) = delete;
    SpatialSearch &operator=(SpatialSearch const &other) = delete;
    SpatialSearch(SpatialSearch &&other) = delete;
    SpatialSearch &operator=(SpatialSearch &&other) = delete;
    SpatialSearch(T const &points, float const cell_size)
    : m_cell_size(cell_size),
      m_points{points}
    {
        update();
    }

    // Note: for optimal results, distance should be lower than
    // the cells' sizes.
    // Note: returned distance is squared.
    std::vector<std::tuple<index_t, index_t, float>> pairs(float const distance) const
    {
        index_t const num_layers = floor(distance / m_cell_size) + 1;
        float const distance2 = distance * distance;
        std::vector<std::tuple<index_t, index_t, float>> pairs_list;

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

    // Note: this function is more efficient than brute-force only if
    // it's used multiple times (e.g. minimum of 30 queries for 30k points).
    // Note: for optimal results, distance should be lower than
    // the cells' sizes.
    std::vector<index_t> query(index_t const query, float const distance) const
    {
        float const distance2 = distance * distance;
        index_t const num_layers = floor(distance / m_cell_size) + 1;
        Point3 point = m_points.col(query);
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

private:
    void update()
    {
        Point3 min_coords = m_points.rowwise().minCoeff();
        Point3 max_coords = m_points.rowwise().maxCoeff();

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

    cell_index_t index(Point3 const &point) const
    {
        return ((point - m_origin) / m_cell_size).array().floor().cast<index_t>();
    }

    cell_index_t clamped_index(Point3 const &point) const
    {
        return index(point).cwiseMax(cell_index_t{1, 1, 1}).cwiseMin(m_max_clamp);
    }

    void find_cells_pairs(std::vector<std::tuple<index_t, index_t, float>> &pairs_list, cell_t const &current, cell_t const &neighbor, float const cutoff2) const
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

            float distance = (m_points.col(i) - m_points.col(j)).squaredNorm();
            if (distance <= cutoff2)
            {
                pairs_list.push_back(std::tuple(i, j, distance));
            }
        }
        }
    }

    float m_cell_size;
    Point3 m_origin;
    cell_index_t m_grid_size;
    cell_index_t m_max_clamp;
    stride_t m_strides;
    T const &m_points;
    std::vector<cell_t> m_cells;
};

} // namespace mol::internal

#endif // SPATIALSEARCH_HPP
