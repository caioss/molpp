#include "matchers.hpp"
#include "auxiliary.hpp"

#include <molpp/Residue.hpp>
#include <molpp/AtomSel.hpp>
#include <molpp/MolError.hpp>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

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
    , name{data.properties().add<Atom, Name>(false)}
    {
        name->value(0) = "A";
        name->value(1) = "B";
        name->value(2) = "C";
    }

    MolData data;
    Aggregate aggregate;
    Aggregate const& const_aggregate;
    Name* name;
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

    EXPECT_FALSE(default_aggregate);
    EXPECT_FALSE(default_const_aggregate);
}

TEST_F(AtomAggregateTest, ValidityOfInvalidIndex)
{
    Aggregate invalid_aggregate(4, 0, &data);
    Aggregate const invalid_const_aggregate(4, 0, &data);

    EXPECT_FALSE(invalid_aggregate);
    EXPECT_FALSE(invalid_const_aggregate);
}

TEST_F(AtomAggregateTest, ValidityOfNullData)
{
    Aggregate null_aggregate(1, 0, nullptr);
    Aggregate const null_const_aggregate(1, 0, nullptr);

    EXPECT_FALSE(null_aggregate);
    EXPECT_FALSE(null_const_aggregate);
}

TEST_F(AtomAggregateTest, Validity)
{
    EXPECT_TRUE(aggregate);
    EXPECT_TRUE(const_aggregate);
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
    EXPECT_THAT(const_aggregate.coords().reshaped(), ElementsAre(1, 1, 1));
    EXPECT_THAT(aggregate.coords().reshaped(), ElementsAre(1, 1, 1));
}

TEST_F(AtomAggregateTest, MutablePositions)
{
    ASSERT_THAT(aggregate.coords().reshaped(), ElementsAre(1, 1, 1));

    aggregate.coords() *= 2;
    EXPECT_THAT(aggregate.coords().reshaped(), ElementsAre(2, 2, 2));
}

TEST_F(AtomAggregateTest, PositionsOnInvalidFrame)
{
    aggregate.set_frame(std::nullopt);

    EXPECT_THROW(aggregate.coords(), MolError);
}

TEST_F(AtomAggregateTest, HasProperty)
{
    ASSERT_THAT(name, NotNull());
    EXPECT_TRUE(aggregate.has<Name>());
    EXPECT_TRUE(const_aggregate.has<Name>());
}

TEST_F(AtomAggregateTest, HasNotProperty)
{
    EXPECT_FALSE(aggregate.has<ResID>());
    EXPECT_FALSE(const_aggregate.has<ResID>());
}

TEST_F(AtomAggregateTest, FetchProperty)
{
    ASSERT_THAT(name, NotNull());
    EXPECT_EQ(aggregate.property<Name>(), name);
    EXPECT_EQ(const_aggregate.property<Name>(), name);
}

TEST_F(AtomAggregateTest, FetchNotPresentProperty)
{
    EXPECT_THAT(aggregate.property<ResID>(), IsNull());
    EXPECT_THAT(const_aggregate.property<ResID>(), IsNull());
}

TEST_F(AtomAggregateTest, GetPropertyValue)
{
    ASSERT_THAT(name, NotNull());

    EXPECT_EQ(aggregate.get<Name>(), "B");
    EXPECT_EQ(const_aggregate.get<Name>(), "B");
}

TEST_F(AtomAggregateTest, GetPropertyValueUsingProperty)
{
    ASSERT_THAT(name, NotNull());

    EXPECT_EQ(aggregate.get(name), "B");
    EXPECT_EQ(const_aggregate.get(name), "B");
}

TEST_F(AtomAggregateTest, SetValueFromGetReference)
{
    ASSERT_THAT(name, NotNull());
    aggregate.get<Name>() = "X";

    EXPECT_EQ(aggregate.get<Name>(), "X");
}

TEST_F(AtomAggregateTest, SetValueFromGetReferenceUsingProperty)
{
    ASSERT_THAT(name, NotNull());
    aggregate.get(name) = "X";

    EXPECT_EQ(aggregate.get(name), "X");
}

TEST_F(AtomAggregateTest, SetPropertyValue)
{
    ASSERT_THAT(name, NotNull());
    aggregate.set<Name>("X");

    EXPECT_EQ(aggregate.get<Name>(), "X");
}

TEST_F(AtomAggregateTest, SetPropertyValueUsingProperty)
{
    ASSERT_THAT(name, NotNull());
    aggregate.set(name, "X");

    EXPECT_EQ(aggregate.get(name), "X");
}
