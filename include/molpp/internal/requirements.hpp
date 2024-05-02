#ifndef MOLPP_INTERNAL_REQUIREMENTS_HPP
#define MOLPP_INTERNAL_REQUIREMENTS_HPP

#include <ranges>
#include <vector>
#include <concepts>

namespace mol::internal
{

template<class T>
concept IndexRange = requires(T t) {
    std::ranges::range<T>&& std::unsigned_integral<typename T::value_type>;
};

} // namespace mol::internal

#endif // MOLPP_INTERNAL_REQUIREMENTS_HPP
