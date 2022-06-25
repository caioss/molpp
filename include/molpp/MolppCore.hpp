#ifndef MOLPPCORE_HPP
#define MOLPPCORE_HPP

#include "MolError.hpp"
#include <Eigen/Dense>
#include <ranges>

using Point3 = Eigen::Vector3f;
using Coord3 = Eigen::Matrix3Xf;

template <class T>
concept SelIndexCompatible = requires(T t)
{
    std::convertible_to<T, std::vector<size_t>>
    && std::ranges::range<T>;
};

#endif // MOLPPCORE_HPP
