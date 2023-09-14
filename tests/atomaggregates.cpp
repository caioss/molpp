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

TEST(AtomAggregates, AtomAggregate) {
    using Aggregate = AtomAggregate<Residue>;
    MolData data = create_moldata(3, 1, 1, 1, 1);

    // Comparison
    EXPECT_TRUE(Aggregate(1, 0, &data) == Aggregate(1, 0, &data));
    EXPECT_FALSE(Aggregate(0, 0, &data) == Aggregate(1, 0, &data));
    EXPECT_FALSE(Aggregate(1, std::nullopt, &data) == Aggregate(1, 0, &data));
    EXPECT_FALSE(Aggregate(1, 0, &data) == Aggregate(1, 0, nullptr));

    Aggregate aggr(1, 0, &data);
    Aggregate const const_aggr(1, 0, &data);

    /*
     * Constructors
     */
    EXPECT_FALSE(Aggregate());
    EXPECT_FALSE(Aggregate(4, 0, &data));
    EXPECT_FALSE(Aggregate(1, 0, nullptr));
    EXPECT_TRUE(aggr);
    EXPECT_TRUE(const_aggr);

    /*
     * Properties
     */
    EXPECT_EQ(aggr.index(), 1);
    EXPECT_EQ(aggr.frame(), 0);
    aggr.set_frame(std::nullopt);
    EXPECT_FALSE(aggr.frame());
    aggr.set_frame(0);
    EXPECT_EQ(aggr.frame(), 0);
    EXPECT_THROW(aggr.set_frame(1), MolError);

    /*
     * Coordinates
     */
    EXPECT_THAT(const_aggr.coords().reshaped(), ElementsAre(1, 1, 1));
    EXPECT_THAT(aggr.coords().reshaped(), ElementsAre(1, 1, 1));
    aggr.coords() *= 2;
    EXPECT_THAT(aggr.coords().reshaped(), ElementsAre(2, 2, 2));

    EXPECT_THROW(Aggregate(1, std::nullopt, &data).coords(), MolError);
}
