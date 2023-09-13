#include "core/PropertyContainer.hpp"
#include <molpp/Atom.hpp>
#include <molpp/Residue.hpp>
#include <molpp/Property.hpp>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace testing;
using namespace mol;
using namespace mol::internal;

class PropertyContainerTest : public ::testing::Test
{
protected:
    PropertyContainerTest()
    {}

    PropertyContainer container;
};

TEST_F(PropertyContainerTest, InitialSize)
{
    container.add<Atom, Mass>(false);
    container.add<Residue, Charge>(false);
    Mass* mass = container.get<Atom, Mass>(0);
    Charge* charge = container.get<Residue, Charge>(0);

    EXPECT_EQ(container.size<Atom>(), 0);
    EXPECT_EQ(container.size<Residue>(), 0);
    EXPECT_EQ(mass->size(), 0);
    EXPECT_EQ(charge->size(), 0);
}

TEST_F(PropertyContainerTest, SetSize)
{
    Mass* mass = container.add<Atom, Mass>(false);
    Charge* charge = container.add<Residue, Charge>(false);

    container.set_size<Atom>(2);
    container.set_size<Residue>(3);

    EXPECT_EQ(container.size<Atom>(), 2);
    EXPECT_EQ(container.size<Residue>(), 3);
    EXPECT_EQ(mass->size(), 2);
    EXPECT_EQ(charge->size(), 3);
}

TEST_F(PropertyContainerTest, GetNotTimeBased)
{
    Mass* from_add = container.add<Atom, Mass>(false);
    Mass* from_get = container.get<Atom, Mass>(0);

    EXPECT_THAT(from_add, NotNull());
    EXPECT_THAT(from_get, NotNull());
    EXPECT_EQ(from_add, from_get);
}

TEST_F(PropertyContainerTest, GetTimeBased)
{
    Mass* from_add = container.add<Atom, Mass>(true);
    Mass* from_get = container.get<Atom, Mass>(0);

    EXPECT_THAT(from_add, IsNull());
    EXPECT_THAT(from_get, IsNull());
}

TEST_F(PropertyContainerTest, GetNotRegistered)
{
    container.add<Atom, Mass>(false);
    Charge* charge = container.get<Atom, Charge>(0);

    EXPECT_THAT(charge, IsNull());
}

TEST_F(PropertyContainerTest, Add)
{
    Mass* time_based = container.add<Atom, Mass>(true);
    Charge* not_time_based = container.add<Atom, Charge>(false);

    EXPECT_THAT(time_based, IsNull());
    EXPECT_THAT(not_time_based, NotNull());
}

TEST_F(PropertyContainerTest, AddTwice)
{
    Mass* first = container.add<Atom, Mass>(false);
    Mass* second = container.add<Atom, Mass>(false);
    Mass* third = container.add<Atom, Mass>(true);

    EXPECT_THAT(first, NotNull());
    EXPECT_THAT(second, IsNull());
    EXPECT_THAT(third, IsNull());
}
