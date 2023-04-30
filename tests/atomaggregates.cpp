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

TEST(AtomAggregates, BaseAtomAggregate) {
    MolData data = create_moldata(3, 1, 1, 1, 0);

    // Comparison
    EXPECT_TRUE(BaseAtomAggregate(1, 0, &data) == BaseAtomAggregate(1, 0, &data));
    EXPECT_FALSE(BaseAtomAggregate(0, 0, &data) == BaseAtomAggregate(1, 0, &data));
    EXPECT_FALSE(BaseAtomAggregate(1, std::nullopt, &data) == BaseAtomAggregate(1, 0, &data));
    EXPECT_FALSE(BaseAtomAggregate(1, 0, &data) == BaseAtomAggregate(1, 0, nullptr));

    /*
     * Properties
     */
    BaseAtomAggregate aggr(1, 0, &data);
    EXPECT_EQ(aggr.index(), 1);
    EXPECT_EQ(aggr.frame(), 0);
}

TEST(AtomAggregates, AtomAggregate) {
    using Aggregate = AtomAggregate<Residue>;
    MolData data = create_moldata(3, 1, 1, 1, 1);

    // Comparison
    EXPECT_TRUE(Aggregate(1, 0, &data) == Aggregate(1, 0, &data));
    EXPECT_FALSE(Aggregate(0, 0, &data) == Aggregate(1, 0, &data));
    EXPECT_FALSE(Aggregate(1, std::nullopt, &data) == Aggregate(1, 0, &data));
    EXPECT_FALSE(Aggregate(1, 0, &data) == Aggregate(1, 0, nullptr));

    /*
     * Properties
     */
    Aggregate aggr(1, 0, &data);
    EXPECT_EQ(aggr.index(), 1);
    EXPECT_EQ(aggr.frame(), 0);

    /*
     * Coordinates
     */
    EXPECT_THAT(aggr.coords().reshaped(), ElementsAre(1, 1, 1));
    aggr.coords() *= 2;
    EXPECT_THAT(aggr.coords().reshaped(), ElementsAre(2, 2, 2));

    EXPECT_THROW(Aggregate(1, std::nullopt, &data).coords(), MolError);
}
