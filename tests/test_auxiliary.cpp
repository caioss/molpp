#include "auxiliary.hpp"
#include <molpp/Atom.hpp>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <span>
#include <initializer_list>

using namespace mol;
using namespace mol::internal;
using namespace testing;

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

class CreateMoldataTest : public testing::Test
{
public:
    CreateMoldataTest()
    : data{create_moldata(3, 2, 2, 1, 2)}
    {
    }

    MolData data;
};

template<IsAtomAggregate Aggregate, IsProperty PropertyType>
class PropertyCreateMoldataTest : public CreateMoldataTest
{
public:
    using value_type = typename PropertyTraits<PropertyType>::value_type;

    PropertyCreateMoldataTest()
    : CreateMoldataTest()
    {
        if constexpr (std::integral<value_type> || std::floating_point<value_type>)
        {
            expected = {0, 1, 2, 3, 4, 5};
        }
        else if constexpr (std::is_same_v<value_type, std::string>)
        {
            expected = {"A", "B", "C", "D", "E", "F"};
        }
        else if constexpr (std::is_same_v<value_type, Point3>)
        {
            expected = {Point3{0, 0, 0}, Point3{1, 1, 1}, Point3{2, 2, 2}, Point3{3, 3, 3}, Point3{4, 4, 4}, Point3{5, 5, 5}};
        }
    }

    void test_property(Frame const frame)
    {
        PropertyType* property = data.properties().get<Aggregate, PropertyType>(frame);
        for (size_t i = 0; i < expected.size(); i++)
        {
            EXPECT_EQ(property->value(i), expected[i]) << "Index " << i;
        }
    }

    std::array<value_type, 6> expected;
};

template<IsProperty PropertyType>
class AtomCreateMoldataTest : public PropertyCreateMoldataTest<Atom, PropertyType>
{
};

TYPED_TEST_SUITE_P(AtomCreateMoldataTest);

TYPED_TEST_P(AtomCreateMoldataTest, Property)
{
    this->test_property(Frame());
}

REGISTER_TYPED_TEST_SUITE_P(AtomCreateMoldataTest, Property);

using AtomProperties = testing::Types<Name, Type, AlternateLocation, InsertionCode, AtomicNumber, Occupancy, TemperatureFactor, Mass, Charge, Radius>;
INSTANTIATE_TYPED_TEST_SUITE_P(Atom, AtomCreateMoldataTest, AtomProperties);

TEST_F(CreateMoldataTest, Bonds)
{
    BondData const& bond_data = data.bonds();
    EXPECT_THAT(bond_data.bonded(0), UnorderedElementsAre(0, 2));
    EXPECT_THAT(bond_data.bonded(1), UnorderedElementsAre());
    EXPECT_THAT(bond_data.bonded(2), UnorderedElementsAre(0, 2, 4));
    EXPECT_THAT(bond_data.bonded(3), UnorderedElementsAre());
    EXPECT_THAT(bond_data.bonded(4), UnorderedElementsAre(2, 4));
    EXPECT_THAT(bond_data.bonded(5), UnorderedElementsAre());
}

TEST_F(CreateMoldataTest, NumberOfFrames)
{
    size_t const num_frames = data.properties().num_frames();
    EXPECT_EQ(num_frames, 2);
}

TEST_F(CreateMoldataTest, Positions)
{
    for (size_t frame = 0; frame < 2; frame++)
    {
        Position* position_property = data.properties().get<Atom, Position>(frame);
        Position::type const& positions = position_property->positions();

        EXPECT_THAT(positions.reshaped(), ElementsAre(0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5));
    }
}
