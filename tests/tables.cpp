#include <molpp/ElementsTable.hpp>
#include "tables/ResiduesTable.hpp"
#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace mol;
using namespace mol::internal;
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

TEST(Tables, Residues) {
    EXPECT_FALSE(RESIDUES_TABLE.contains("UNX"));
    EXPECT_TRUE(RESIDUES_TABLE.contains("PHE"));
    EXPECT_THAT(RESIDUES_TABLE.max_atoms(), Gt(0));

    EXPECT_EQ(RESIDUES_TABLE["GLY"].atoms.size(), 10);
    auto const &glycine = RESIDUES_TABLE["GLY"];
    EXPECT_EQ(glycine.atoms.size(), 10);
    EXPECT_EQ(glycine.atom_index("H2"), 6);
    EXPECT_EQ(glycine.atom_index("UNK"), -1);
    EXPECT_EQ(glycine.bonds.size(), 9);

    auto const &bond = glycine.bonds[0];
    EXPECT_FALSE(bond.aromatic);
    EXPECT_EQ(bond.order, 1);
    EXPECT_EQ(bond.atom1, 0);
    EXPECT_EQ(bond.atom2, 1);
}
