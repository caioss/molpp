#include <molpp/ElementsTable.hpp>
#include "tables/ResiduesTable.hpp"
#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace mol;
using namespace mol::internal;
using namespace testing;

TEST(Tables, Elements) {
    ElementsTable const& table = ELEMENTS_TABLE();

    for (int atomic = 0; atomic < 119; ++atomic)
    {
        EXPECT_EQ(table.atomic_number(atomic), atomic) << atomic;
    }

    EXPECT_FLOAT_EQ(table.mass(0), 0);
    EXPECT_FLOAT_EQ(table.mass(6), 12.011);
    EXPECT_FLOAT_EQ(table.mass(118), 294);

    EXPECT_FLOAT_EQ(table.covalent_radius(0), 0);
    EXPECT_FLOAT_EQ(table.covalent_radius(6), 0.77);
    EXPECT_FLOAT_EQ(table.covalent_radius(118), 0);
    EXPECT_FLOAT_EQ(table.max_covalent_radius(), 2.25);

    EXPECT_FLOAT_EQ(table.VDW_radius(0), 0);
    EXPECT_FLOAT_EQ(table.VDW_radius(6), 1.7);
    EXPECT_FLOAT_EQ(table.VDW_radius(118), 0);
    EXPECT_FLOAT_EQ(table.max_VDW_radius(), 3.0);

    EXPECT_EQ(table.symbol(0), "Xx");
    EXPECT_EQ(table.symbol(6), "C");
    EXPECT_EQ(table.symbol(118), "Uuo");

    EXPECT_EQ(table.name(0), "Dummy");
    EXPECT_EQ(table.name(6), "Carbon");
    EXPECT_EQ(table.name(118), "Ununoctium");

    auto carbon = table(6);
    EXPECT_EQ(carbon.atomic_number, 6);
    EXPECT_FLOAT_EQ(carbon.mass, 12.011);
    EXPECT_FLOAT_EQ(carbon.covalent_radius, 0.77);
    EXPECT_FLOAT_EQ(carbon.VDW_radius, 1.7);
    EXPECT_EQ(carbon.symbol, "C");
    EXPECT_EQ(carbon.name, "Carbon");

}

TEST(Tables, Residues) {
    ResiduesTable const& table = RESIDUES_TABLE();

    EXPECT_FALSE(table.contains("UNX"));
    EXPECT_TRUE(table.contains("PHE"));
    EXPECT_THAT(table.max_atoms(), Gt(0));

    EXPECT_EQ(table["GLY"].atoms.size(), 10);
    auto const &glycine = table["GLY"];
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
