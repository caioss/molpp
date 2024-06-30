#include "matchers.hpp"
#include "utils.hpp"

#include <molpp/Residue.hpp>
#include <molpp/AtomSel.hpp>
#include <molpp/Error.hpp>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <ranges>

using namespace mol;
using namespace mol::internal;
using namespace testing;

class EntityTest : public ::testing::Test
{
public:
    EntityTest()
    : data{create_moldata(3, 1, 1, 1, 1)}
    , entity{1, 0, data}
    , const_entity{entity}
    {
    }

    MolData data;
    Entity entity;
    Entity const& const_entity;
};

TEST_F(EntityTest, EqualityOperator)
{
    EXPECT_TRUE(Entity(1, 0, data) == Entity(1, 0, data));
    EXPECT_FALSE(Entity(0, 0, data) == Entity(1, 0, data));
    EXPECT_FALSE(Entity(1, std::nullopt, data) == Entity(1, 0, data));
}

TEST_F(EntityTest, Index)
{
    EXPECT_EQ(entity.index(), 1);
    EXPECT_EQ(const_entity.index(), 1);
}

TEST_F(EntityTest, Indices)
{
    EXPECT_THAT(entity.indices(), ElementsAre(1));
    EXPECT_THAT(const_entity.indices(), ElementsAre(1));
}

TEST_F(EntityTest, Frame)
{
    EXPECT_EQ(entity.frame(), 0);
    EXPECT_EQ(const_entity.frame(), 0);
}

TEST_F(EntityTest, SetNullFrame)
{
    entity.set_frame(std::nullopt);

    EXPECT_FALSE(entity.frame());
}

TEST_F(EntityTest, ChangeFrame)
{
    entity.set_frame(std::nullopt);
    entity.set_frame(0);

    EXPECT_EQ(entity.frame(), 0);
}

TEST_F(EntityTest, SetOutOfRangeFrame)
{
    EXPECT_THROW(entity.set_frame(1), Error);
}
