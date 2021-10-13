#ifndef SPATIALSEARCH_HPP
#define SPATIALSEARCH_HPP

#include <molpp/MolppCore.hpp>
#include <vector>
#include <utility>
#include <cstdlib>

namespace mol::internal {

class SpatialSearch {
public:
    using index_t = ptrdiff_t;

private:
    using stride_t = Eigen::RowVector3<index_t>;
    using cell_index_t = Eigen::Vector3<index_t>;
    using point_t = Eigen::Vector3f;
    using point_data_t = Eigen::Matrix3Xf;
    using cell_t = std::vector<index_t>;

public:
    SpatialSearch() = delete;
    SpatialSearch(SpatialSearch const &other) = delete;
    SpatialSearch &operator=(SpatialSearch const &other) = delete;
    SpatialSearch(SpatialSearch &&other) = delete;
    SpatialSearch &operator=(SpatialSearch &&other) = delete;
    SpatialSearch(Eigen::Ref<point_data_t> const &points, float const cell_size);

    // Note: for optimal results, distance should be lower than
    // the cells' sizes.
    std::vector<std::pair<index_t, index_t>> pairs(float const distance) const;

    // Note: this function is more efficient than brute-force only if
    // it's used multiple times (e.g. minimum of 30 queries for 30k points).
    // Note: for optimal results, distance should be lower than
    // the cells' sizes.
    std::vector<index_t> query(index_t const query, float const distance) const;

private:
    void update();
    cell_index_t index(Eigen::Ref<point_t> const &point) const;
    cell_index_t clamped_index(Eigen::Ref<point_t> const &point) const;
    void find_cells_pairs(std::vector<std::pair<index_t, index_t>> &pairs_list, cell_t const &current, cell_t const &neighbor, float const distance2) const;

    float m_cell_size;
    point_t m_origin;
    cell_index_t m_grid_size;
    cell_index_t m_max_clamp;
    stride_t m_strides;
    Eigen::Ref<point_data_t> m_points;
    std::vector<cell_t> m_cells;
};

} // namespace mol::internal

#endif // SPATIALSEARCH_HPP
