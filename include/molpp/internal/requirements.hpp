#ifndef REQUIREMENTS_HPP
#define REQUIREMENTS_HPP

#include <ranges>
#include <vector>
#include <concepts>

namespace mol::internal
{

template <class T>
concept IndexRange = requires(T t)
{
    std::ranges::range<T>
    && std::unsigned_integral<typename T::value_type>;
};

template <class Derived>
class AtomAggregate;

template <class T>
concept IsAtomAggregate = requires(T t)
{
    std::derived_from<T, AtomAggregate<T>>;
    {t.atom_indices()} -> std::convertible_to<std::vector<index_t>>;
};

} // namespace mol::internal

#endif // REQUIREMENTS_HPP
