#ifndef MATH_HPP
#define MATH_HPP

#include <concepts>

namespace mol::internal {

// Floating point numbers comparisons.
// Adapted from: The Art of Computer Programming, Donald KnuthV
template <std::floating_point T>
bool approximately_equal(T a, T b, T epsilon)
{
    return std::abs(a - b) <= (std::max(std::abs(a), std::abs(b)) * epsilon);
}

template <std::floating_point T>
bool essentially_equal(T a, T b, T epsilon)
{
    return std::abs(a - b) <= (std::min(std::abs(a), std::abs(b)) * epsilon);
}

} // namespace mol::internal

#endif // MATH_HPP
