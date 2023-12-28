#include <molpp/Property.hpp>
#include <molpp/MolppCore.hpp>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <array>

using namespace testing;
using namespace mol;

template <IsProperty Type>
struct PropertyTraits
{
    using property_type = Type;
    using value_type = Type::value_type;
};

template <>
struct PropertyTraits<Position>
{
    using property_type = Position;
    using value_type = Point3;
};

template<IsProperty Type>
class PropertyTestValues
{
public:
    using value_type = typename PropertyTraits<Type>::value_type;

    PropertyTestValues()
    {
        if constexpr (std::integral<value_type> || std::floating_point<value_type>)
        {
            expected = {1, 2, 3, 4, 5, 6};
        }
        else if constexpr (std::is_same_v<value_type, std::string>)
        {
            expected = {"a", "b", "c", "d", "e", "f"};
        }
        else if constexpr (std::is_same_v<value_type, Point3>)
        {
            expected = {Point3{0, 0, 0}, Point3{1, 1, 1}, Point3{2, 2, 2}, Point3{3, 3, 3}, Point3{4, 4, 4}, Point3{5, 5, 5}};
        }
    }

    std::array<value_type, 6> expected;
};

template<IsProperty Type>
class PropertyTest : public Type, public PropertyTestValues<Type>, public testing::Test
{
};

TYPED_TEST_SUITE_P(PropertyTest);

TYPED_TEST_P(PropertyTest, ResizeOnce)
{
    this->resize(2);

    EXPECT_EQ(this->size(), 2);
}

TYPED_TEST_P(PropertyTest, ResizeTwice)
{
    this->resize(2);
    this->resize(4);

    EXPECT_EQ(this->size(), 4);
}

TYPED_TEST_P(PropertyTest, SetValues)
{
    this->resize(2);
    this->value(0) = this->expected[0];
    this->value(1) = this->expected[1];

    EXPECT_EQ(this->size(), 2);
    EXPECT_EQ(this->value(0), this->expected[0]);
    EXPECT_EQ(this->value(1), this->expected[1]);
}

TYPED_TEST_P(PropertyTest, ResizeTwiceAndSetValues)
{
    this->resize(4);
    this->resize(2);
    this->value(0) = this->expected[0];
    this->value(1) = this->expected[1];

    EXPECT_EQ(this->size(), 2);
    EXPECT_EQ(this->value(0), this->expected[0]);
    EXPECT_EQ(this->value(1), this->expected[1]);
}

REGISTER_TYPED_TEST_SUITE_P(PropertyTest, ResizeOnce, ResizeTwice, SetValues, ResizeTwiceAndSetValues);

using Properties = testing::Types<Occupancy, TemperatureFactor, Mass, Charge, Radius, AtomicNumber, ResID, Name, Type, AlternateLocation, InsertionCode, ResName, Position>;
INSTANTIATE_TYPED_TEST_SUITE_P(VectorProperties, PropertyTest, Properties);

using PositionPropertyTest = PropertyTest<Position>;

TEST_F(PositionPropertyTest, Positions)
{
    this->resize(this->expected.size());
    for (size_t i = 0; i < this->expected.size(); i++)
    {
        this->value(i) = this->expected[i];
    }

    for (size_t i = 0; i < this->expected.size(); i++)
    {
        EXPECT_EQ(this->value(i), this->expected[i]) << i;
    }
}
