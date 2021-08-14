#include "ElementsTable.hpp"
#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace mol;
using namespace testing;

TEST(Tables, Elements) {
    for (int atomic = 0; atomic < 119; ++atomic)
    {
        EXPECT_EQ(ELEMENTS_TABLE.atomic_number(atomic), atomic) << atomic;
    }

    EXPECT_FLOAT_EQ(ELEMENTS_TABLE.mass(0), 0);
    EXPECT_FLOAT_EQ(ELEMENTS_TABLE.mass(6), 12.011);
    EXPECT_FLOAT_EQ(ELEMENTS_TABLE.mass(118), 294);

    EXPECT_FLOAT_EQ(ELEMENTS_TABLE.covalent_radius(0), 0);
    EXPECT_FLOAT_EQ(ELEMENTS_TABLE.covalent_radius(6), 0.77);
    EXPECT_FLOAT_EQ(ELEMENTS_TABLE.covalent_radius(118), 0);

    EXPECT_FLOAT_EQ(ELEMENTS_TABLE.VDW_radius(0), 0);
    EXPECT_FLOAT_EQ(ELEMENTS_TABLE.VDW_radius(6), 1.7);
    EXPECT_FLOAT_EQ(ELEMENTS_TABLE.VDW_radius(118), 0);

    EXPECT_EQ(ELEMENTS_TABLE.symbol(0), "Xx");
    EXPECT_EQ(ELEMENTS_TABLE.symbol(6), "C");
    EXPECT_EQ(ELEMENTS_TABLE.symbol(118), "Uuo");

    EXPECT_EQ(ELEMENTS_TABLE.name(0), "Dummy");
    EXPECT_EQ(ELEMENTS_TABLE.name(6), "Carbon");
    EXPECT_EQ(ELEMENTS_TABLE.name(118), "Ununoctium");

    auto carbon = ELEMENTS_TABLE(6);
    EXPECT_EQ(carbon.atomic_number, 6);
    EXPECT_FLOAT_EQ(carbon.mass, 12.011);
    EXPECT_FLOAT_EQ(carbon.covalent_radius, 0.77);
    EXPECT_FLOAT_EQ(carbon.VDW_radius, 1.7);
    EXPECT_EQ(carbon.symbol, "C");
    EXPECT_EQ(carbon.name, "Carbon");
}
