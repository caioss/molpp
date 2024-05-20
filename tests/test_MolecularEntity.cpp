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

class AtomAggregateTest : public ::testing::Test
{
public:
    AtomAggregateTest()
    : data{create_moldata(3, 1, 1, 1, 1)}
    , aggregate{1, 0, data}
    , const_aggregate{aggregate}
    {
    }

    MolData data;
    MolecularEntity aggregate;
    MolecularEntity const& const_aggregate;
};

TEST_F(AtomAggregateTest, EqualityOperator)
{
    EXPECT_TRUE(MolecularEntity(1, 0, data) == MolecularEntity(1, 0, data));
    EXPECT_FALSE(MolecularEntity(0, 0, data) == MolecularEntity(1, 0, data));
    EXPECT_FALSE(MolecularEntity(1, std::nullopt, data) == MolecularEntity(1, 0, data));
}

TEST_F(AtomAggregateTest, Index)
{
    EXPECT_EQ(aggregate.index(), 1);
    EXPECT_EQ(const_aggregate.index(), 1);
}

TEST_F(AtomAggregateTest, Frame)
{
    EXPECT_EQ(aggregate.frame(), 0);
    EXPECT_EQ(const_aggregate.frame(), 0);
}

TEST_F(AtomAggregateTest, SetNullFrame)
{
    aggregate.set_frame(std::nullopt);

    EXPECT_FALSE(aggregate.frame());
}

TEST_F(AtomAggregateTest, ChangeFrame)
{
    aggregate.set_frame(std::nullopt);
    aggregate.set_frame(0);

    EXPECT_EQ(aggregate.frame(), 0);
}

TEST_F(AtomAggregateTest, SetOutOfRangeFrame)
{
    EXPECT_THROW(aggregate.set_frame(1), MolError);
}
