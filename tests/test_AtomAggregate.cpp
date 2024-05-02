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
    , aggregate{1, 0, &data}
    , const_aggregate{aggregate}
    {
    }

    MolData data;
    AtomAggregate aggregate;
    AtomAggregate const& const_aggregate;
};

TEST_F(AtomAggregateTest, EqualityOperator)
{
    EXPECT_TRUE(AtomAggregate(1, 0, &data) == AtomAggregate(1, 0, &data));
    EXPECT_FALSE(AtomAggregate(0, 0, &data) == AtomAggregate(1, 0, &data));
    EXPECT_FALSE(AtomAggregate(1, std::nullopt, &data) == AtomAggregate(1, 0, &data));
    EXPECT_FALSE(AtomAggregate(1, 0, &data) == AtomAggregate(1, 0, nullptr));
}

TEST_F(AtomAggregateTest, ValidityOfDefaultConstructed)
{
    AtomAggregate default_aggregate{};
    AtomAggregate const default_const_aggregate{};

    EXPECT_FALSE(default_aggregate.is_valid());
    EXPECT_FALSE(default_const_aggregate.is_valid());
}

TEST_F(AtomAggregateTest, ValidityOfNullData)
{
    AtomAggregate null_aggregate(1, 0, nullptr);
    AtomAggregate const null_const_aggregate(1, 0, nullptr);

    EXPECT_FALSE(null_aggregate.is_valid());
    EXPECT_FALSE(null_const_aggregate.is_valid());
}

TEST_F(AtomAggregateTest, Validity)
{
    EXPECT_TRUE(aggregate.is_valid());
    EXPECT_TRUE(const_aggregate.is_valid());
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
