#include <molpp/Property.hpp>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <array>

using namespace testing;
using namespace mol;

template<IsProperty T>
class PropertyTest : public testing::Test
{
protected:
    PropertyTest()
    : property()
    {
        if constexpr (std::integral<typename T::value_type> || std::floating_point<typename T::value_type>)
        {
            values = {1, 2, 3, 4, 5, 6};
        }
        else if constexpr (std::is_same_v<typename T::value_type, std::string>)
        {
            values = {"a", "b", "c", "d", "e", "f"};
        }
    }

    T property;
    std::array<typename T::value_type, 6> values;
};

TYPED_TEST_SUITE_P(PropertyTest);

TYPED_TEST_P(PropertyTest, ResizeOnce)
{
    this->property.resize(2);

    EXPECT_EQ(this->property.size(), 2);
}

TYPED_TEST_P(PropertyTest, ResizeTwice)
{
    this->property.resize(2);
    this->property.resize(4);

    EXPECT_EQ(this->property.size(), 4);
}

TYPED_TEST_P(PropertyTest, SetValues)
{
    this->property.resize(2);
    this->property.value(0) = this->values[0];
    this->property.value(1) = this->values[1];

    EXPECT_EQ(this->property.size(), 2);
    EXPECT_EQ(this->property.value(0), this->values[0]);
    EXPECT_EQ(this->property.value(1), this->values[1]);
}

TYPED_TEST_P(PropertyTest, ResizeTwiceWithValues)
{
    this->property.resize(4);
    this->property.value(0) = this->values[0];
    this->property.value(1) = this->values[1];
    this->property.value(2) = this->values[2];
    this->property.value(3) = this->values[3];
    this->property.resize(2);

    EXPECT_EQ(this->property.size(), 2);
    EXPECT_EQ(this->property.value(0), this->values[0]);
    EXPECT_EQ(this->property.value(1), this->values[1]);
}

REGISTER_TYPED_TEST_SUITE_P(PropertyTest, ResizeOnce, ResizeTwice, SetValues, ResizeTwiceWithValues);

using ContainerPropertyTypes = testing::Types<ContainerProperty<int>, ContainerProperty<std::string>>;
INSTANTIATE_TYPED_TEST_SUITE_P(ContainerProperties, PropertyTest, ContainerPropertyTypes);

using PropertyTypes = testing::Types<Occupancy, TemperatureFactor, Mass, Charge, Radius, AtomicNumber, ResID, Name, Type, AlternateLocation, InsertionCode, ResName>;
INSTANTIATE_TYPED_TEST_SUITE_P(VectorProperties, PropertyTest, PropertyTypes);
