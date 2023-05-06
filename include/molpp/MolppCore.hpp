#ifndef MOLPPCORE_HPP
#define MOLPPCORE_HPP

#include "MolError.hpp"
#include <Eigen/Dense>
#include <optional>

namespace mol {

// Substructure indexing (Atom, Residue, etc...)
using index_t = size_t;
using position_t = float;

using Point3 = Eigen::Vector<position_t, 3>;
using Coord3 = Eigen::Matrix<position_t, 3, Eigen::Dynamic>;
using Coord2 = Eigen::Matrix<position_t, 2, Eigen::Dynamic>;

using Frame = std::optional<size_t>;

}; // namespace mol

#endif // MOLPPCORE_HPP
