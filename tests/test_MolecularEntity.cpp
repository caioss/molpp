#include "matchers.hpp"
#include "auxiliary.hpp"

#include <molpp/Residue.hpp>
#include <molpp/AtomSel.hpp>
#include <molpp/MolError.hpp>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <ranges>

using namespace mol;
using namespace mol::internal;
using namespace testing;

class MolecularEntityTest : public ::testing::Test
{
public:
    MolecularEntityTest()
    : data{create_moldata(3, 1, 1, 1, 1)}
    , entity{1, 0, data}
    , const_entity{entity}
    {
    }

    MolData data;
    MolecularEntity entity;
    MolecularEntity const& const_entity;
};

TEST_F(MolecularEntityTest, EqualityOperator)
{
    EXPECT_TRUE(MolecularEntity(1, 0, data) == MolecularEntity(1, 0, data));
    EXPECT_FALSE(MolecularEntity(0, 0, data) == MolecularEntity(1, 0, data));
    EXPECT_FALSE(MolecularEntity(1, std::nullopt, data) == MolecularEntity(1, 0, data));
}

TEST_F(MolecularEntityTest, Index)
{
    EXPECT_EQ(entity.index(), 1);
    EXPECT_EQ(const_entity.index(), 1);
}

TEST_F(MolecularEntityTest, Indices)
{
    EXPECT_THAT(entity.indices(), ElementsAre(1));
    EXPECT_THAT(const_entity.indices(), ElementsAre(1));
}

TEST_F(MolecularEntityTest, Frame)
{
    EXPECT_EQ(entity.frame(), 0);
    EXPECT_EQ(const_entity.frame(), 0);
}

TEST_F(MolecularEntityTest, SetNullFrame)
{
    entity.set_frame(std::nullopt);

    EXPECT_FALSE(entity.frame());
}

TEST_F(MolecularEntityTest, ChangeFrame)
{
    entity.set_frame(std::nullopt);
    entity.set_frame(0);

    EXPECT_EQ(entity.frame(), 0);
}

TEST_F(MolecularEntityTest, SetOutOfRangeFrame)
{
    EXPECT_THROW(entity.set_frame(1), MolError);
}
