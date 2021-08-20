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
    auto data = aux_moldata();

    // Comparison
    EXPECT_TRUE(BaseAtomAggregate(1, 0, data) == BaseAtomAggregate(1, 0, data));
    EXPECT_FALSE(BaseAtomAggregate(0, 0, data) == BaseAtomAggregate(1, 0, data));
    EXPECT_FALSE(BaseAtomAggregate(1, std::nullopt, data) == BaseAtomAggregate(1, 0, data));
    EXPECT_FALSE(BaseAtomAggregate(1, 0, data) == BaseAtomAggregate(1, 0, nullptr));

    /*
     * Properties
     */
    BaseAtomAggregate aggr(1, 0, data);
    EXPECT_EQ(aggr.index(), 1);
    EXPECT_EQ(aggr.frame(), 0);
}

TEST(AtomAggregates, AtomAggregate) {
    using Aggregate = AtomAggregate<Residue>;
    auto data = aux_moldata();

    // Comparison
    EXPECT_TRUE(Aggregate(1, 0, data) == Aggregate(1, 0, data));
    EXPECT_FALSE(Aggregate(0, 0, data) == Aggregate(1, 0, data));
    EXPECT_FALSE(Aggregate(1, std::nullopt, data) == Aggregate(1, 0, data));
    EXPECT_FALSE(Aggregate(1, 0, data) == Aggregate(1, 0, nullptr));

    /*
     * Properties
     */
    Aggregate aggr(1, 0, data);
    EXPECT_EQ(aggr.index(), 1);
    EXPECT_EQ(aggr.frame(), 0);

    /*
     * Coordinates
     */
    EXPECT_THAT(aggr.coords().reshaped(), ElementsAre(2, 5, 8));
    aggr.coords() *= 2;
    EXPECT_THAT(aggr.coords().reshaped(), ElementsAre(4, 10, 16));

    EXPECT_THROW(Aggregate(1, std::nullopt, data).coords(), MolError);

    /*
     * Conversions
     */
    auto sel = aggr.atoms();
    ASSERT_THAT(sel, NotNull());
    ASSERT_EQ(sel->size(), 1);
    EXPECT_EQ(sel->frame(), aggr.frame());
    EXPECT_EQ((*sel)[0].index(), 1);
    EXPECT_EQ((*sel)[0].residue_id(), 1);
}
