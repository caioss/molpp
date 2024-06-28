#ifndef AUXILIARY_HPP
#define AUXILIARY_HPP

#include <molpp/internal/MolData.hpp>

#include <vector>
#include <type_traits>

mol::internal::MolData create_moldata(size_t const num_res, size_t const num_res_atoms, size_t const num_chains, size_t const num_segments, size_t const num_frames);

template<std::ranges::range Range>
constexpr auto view2vector(Range&& r)
{
    using elem_t = std::decay_t<std::ranges::range_value_t<Range>>;
    return std::vector<elem_t>{r.begin(), r.end()};
}

namespace
{

//! Helper to get information about member functions
template<class Return = void, class Class = void, class... Args>
struct MemberFunctionTraitsDetails
{
    using class_type = Class;
    using class_base_type = std::remove_cvref_t<Class>;
    using return_type = Return;
    using return_base_type = std::remove_cvref_t<Return>;
    using arguments_type = std::tuple<Args...>;

    explicit MemberFunctionTraitsDetails(auto)
    {}
};

template<class Return, class Class, class... Args>
MemberFunctionTraitsDetails(Return (Class::*)(Args...)) -> MemberFunctionTraitsDetails<Return, Class, Args...>;

template<class Return, class Class, class... Args>
MemberFunctionTraitsDetails(Return (Class::*)(Args...) const) -> MemberFunctionTraitsDetails<Return, Class, Args...>;

} // namespace

//! Traits for member functions
template<class MFP>
using MemberFunctionTraits = decltype(MemberFunctionTraitsDetails((MFP){}));

#endif // AUXILIARY_HPP
