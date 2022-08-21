#ifndef MOLPPCORE_HPP
#define MOLPPCORE_HPP

#include "MolError.hpp"
#include <Eigen/Dense>
#include <ranges>

namespace mol {

using position_t = float;

using Point3 = Eigen::Vector<position_t, 3>;
using Coord3 = Eigen::Matrix<position_t, 3, Eigen::Dynamic>;
using Coord2 = Eigen::Matrix<position_t, 2, Eigen::Dynamic>;

template <class T>
concept SelIndexCompatible = requires(T t)
{
    std::convertible_to<T, std::vector<size_t>>
    && std::ranges::range<T>;
};

}; // namespace mol

#endif // MOLPPCORE_HPP
