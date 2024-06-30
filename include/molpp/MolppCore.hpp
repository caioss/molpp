#ifndef MOLPP_MOLPPCORE_HPP
#define MOLPP_MOLPPCORE_HPP

#include "MolError.hpp"

#include <Eigen/Dense>

#include <ranges>
#include <optional>
#include <concepts>

namespace mol
{

// Substructure indexing (Atom, Residue, etc...)
using index_t = size_t;
using position_t = float;

enum class MolecularEntityCategory
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

using Point3 = Eigen::Vector<position_t, 3>;
using Coord3 = Eigen::Matrix<position_t, 3, Eigen::Dynamic>;
using Coord2 = Eigen::Matrix<position_t, 2, Eigen::Dynamic>;

using Frame = std::optional<size_t>;

namespace internal
{

template<class Range>
concept IndexRange = std::ranges::range<Range> && std::same_as<std::ranges::range_value_t<Range>, index_t>;

} // namespace internal

}; // namespace mol

#endif // MOLPP_MOLPPCORE_HPP
