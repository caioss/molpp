#ifndef REQUIREMENTS_HPP
#define REQUIREMENTS_HPP

#include <ranges>
#include <concepts>

namespace mol::internal
{

template <class T>
concept IndexRange = requires(T t)
{
    std::ranges::range<T>
    && std::unsigned_integral<typename T::value_type>;
};

} // namespace mol::internal

#endif // REQUIREMENTS_HPP
