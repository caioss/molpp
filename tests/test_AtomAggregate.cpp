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
    using Aggregate = AtomAggregate<Atom>;

    AtomAggregateTest()
    : data{create_moldata(3, 1, 1, 1, 1)}
    , aggregate{1, 0, &data}
    , const_aggregate{aggregate}
    {
    }

    MolData data;
    Aggregate aggregate;
    Aggregate const& const_aggregate;
};

TEST_F(AtomAggregateTest, EqualityOperator)
{
    EXPECT_TRUE(Aggregate(1, 0, &data) == Aggregate(1, 0, &data));
    EXPECT_FALSE(Aggregate(0, 0, &data) == Aggregate(1, 0, &data));
    EXPECT_FALSE(Aggregate(1, std::nullopt, &data) == Aggregate(1, 0, &data));
    EXPECT_FALSE(Aggregate(1, 0, &data) == Aggregate(1, 0, nullptr));
}

TEST_F(AtomAggregateTest, ValidityOfDefaultConstructed)
{
    Aggregate default_aggregate{};
    Aggregate const default_const_aggregate{};

    EXPECT_FALSE(default_aggregate.is_valid());
    EXPECT_FALSE(default_const_aggregate.is_valid());
}

TEST_F(AtomAggregateTest, ValidityOfNullData)
{
    Aggregate null_aggregate(1, 0, nullptr);
    Aggregate const null_const_aggregate(1, 0, nullptr);

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

TEST_F(AtomAggregateTest, Positions)
{
    EXPECT_THAT(aggregate.coords().reshaped(), ElementsAre(1, 1, 1));
    EXPECT_THAT(const_aggregate.coords().reshaped(), ElementsAre(1, 1, 1));
}

TEST_F(AtomAggregateTest, Bonds)
{
    std::vector<std::pair<index_t, index_t>> bonds_indices;
    for (std::shared_ptr<Bond> bond : aggregate.bonds())
    {
        bonds_indices.push_back(std::make_pair<index_t, index_t>(bond->atom1(), bond->atom2()));
    }

    EXPECT_THAT(bonds_indices, UnorderedElementsAre(Pair(0, 1), Pair(1, 2)));
}
