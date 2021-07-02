#include "files.hpp"
#include "matchers.hpp"
#include "core/AtomData.hpp"
#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace mol;
using namespace mol::internal;
using namespace testing;

#if 0
TEST_F(MolFiles, bond_graph) {
}
#endif

TEST(Bonds, bond_data) {
    BondData bond;

    EXPECT_EQ(bond.order(), BondData::Unknown);
    bond.set_order(BondData::Single);
    EXPECT_EQ(bond.order(), BondData::Single);

    EXPECT_TRUE(bond.guessed());
    bond.set_guessed(false);
    EXPECT_FALSE(bond.guessed());

    EXPECT_TRUE(bond.guessed_order());
    bond.set_guessed_order(false);
    EXPECT_FALSE(bond.guessed_order());
}

#if 0
TEST_F(MolFiles, bond_atom) {
}

TEST_F(MolFiles, bond_atom_sel) {
}

TEST_F(MolFiles, bond) {
}
#endif
