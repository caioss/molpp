#ifndef MOLPP_COMMON_HPP
#define MOLPP_COMMON_HPP

#include "Error.hpp"

#include <Eigen/Dense>

#include <ranges>
#include <optional>
#include <concepts>

namespace mol
{

//! Substructure indexing (Atom, Residue, etc...).
using index_t = size_t;

//! Type of positions.
using position_t = float;

//! Entity category tag.
enum class EntityCategory
{
    Atom,
    Residue,
    Chain,
    Segment,
    Custom0,
    Custom1,
    Custom2,
    Custom3,
    Custom4
};

//! Storage of a 3D point.
using Point3 = Eigen::Vector<position_t, 3>;

//! Storage of a set of 3D points.
using Positions3 = Eigen::Matrix<position_t, 3, Eigen::Dynamic>;

//! Storage of a set of 2D points.
using Positions2 = Eigen::Matrix<position_t, 2, Eigen::Dynamic>;

//! Frame number that may not be defined.
using Frame = std::optional<size_t>;

namespace internal
{

//! Concept for a range of indices.
template<class Range>
concept IndexRange = std::ranges::range<Range> && std::same_as<std::ranges::range_value_t<Range>, index_t>;

} // namespace internal

}; // namespace mol

#endif // MOLPP_COMMON_HPP
