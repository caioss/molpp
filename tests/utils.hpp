#ifndef TESTS_UTILS_HPP
#define TESTS_UTILS_HPP

#include <molpp/internal/MolData.hpp>

#include <gtest/gtest.h>

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

//! Test an entity's property against a set of expected values
template<class Property, class ArrayType = MemberFunctionTraits<Property>::return_base_type, size_t Size>
void test_property(Property property, std::array<ArrayType, Size> const expected, mol::internal::MolData& data)
{
    for (size_t i = 0; i < Size; i++)
    {
        using PropertyTraits = MemberFunctionTraits<Property>;
        using Entity = PropertyTraits::class_base_type;

        Entity entity(i, std::nullopt, data);
        Entity const const_entity(i, std::nullopt, data);

        typename PropertyTraits::return_type const value = std::invoke(property, &entity);
        typename PropertyTraits::return_type const const_value = std::invoke(property, &const_entity);

        if constexpr (std::is_same_v<typename PropertyTraits::return_base_type, float>)
        {
            EXPECT_FLOAT_EQ(value, expected[i]) << i << " (non-const)";
            EXPECT_FLOAT_EQ(const_value, expected[i]) << i << " (const)";
        }
        else if constexpr (std::is_same_v<typename PropertyTraits::return_base_type, double>)
        {
            EXPECT_DOUBLE_EQ(value, expected[i]) << i << " (non-const)";
            EXPECT_DOUBLE_EQ(const_value, expected[i]) << i << " (const)";
        }
        else
        {
            EXPECT_EQ(value, expected[i]) << i << " (non-const)";
            EXPECT_EQ(const_value, expected[i]) << i << " (const)";
        }
    }
}

#endif // TESTS_UTILS_HPP
