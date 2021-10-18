#include "files.hpp"
#include "matchers.hpp"
#include "core/MolData.hpp"
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <vector>

using namespace mol;
using namespace mol::internal;
using namespace testing;

TEST(Bonds, BondData) {
    BondData bond_graph(5);

    EXPECT_TRUE(bond_graph.incomplete());
    bond_graph.set_incomplete(false);
    EXPECT_FALSE(bond_graph.incomplete());

    EXPECT_THAT(bond_graph.add_bond(1, 2), NotNull());
    EXPECT_EQ(bond_graph.add_bond(2, 1), bond_graph.bond(1, 2));
    EXPECT_THAT(bond_graph.add_bond(3, 1), NotNull());
    EXPECT_THAT(bond_graph.add_bond(4, 3), NotNull());
    EXPECT_THAT(bond_graph.add_bond(1, 1), IsNull());

    EXPECT_THAT(bond_graph.bonded(0), UnorderedElementsAre());
    EXPECT_THAT(bond_graph.bonded(1), UnorderedElementsAre(1, 2, 3));
    EXPECT_THAT(bond_graph.bonded(2), UnorderedElementsAre(1, 2));
    EXPECT_THAT(bond_graph.bonded(3), UnorderedElementsAre(1, 3, 4));
    EXPECT_THAT(bond_graph.bonded(4), ElementsAre(3, 4));

    size_t indices[3] = {1, 0, 4};
    auto bonded = bond_graph.bonds(indices, indices + 3);
    EXPECT_EQ(bonded.size(), 3);
    std::vector<size_t> bond_indices;
    bond_indices.reserve(6);
    for (auto b : bonded)
    {
        bond_indices.push_back(b->atom1());
        bond_indices.push_back(b->atom2());
    }
    EXPECT_THAT(bond_indices, UnorderedElementsAre(1, 2, 3, 1, 3, 4));
    EXPECT_THAT(bond_graph.bonded(indices, indices + 3), UnorderedElementsAre(1, 2, 3, 4));

    bonded = bond_graph.bonds(1);
    EXPECT_EQ(bonded.size(), 2);
    bond_indices.clear();
    bond_indices.reserve(4);
    for (auto b : bonded)
    {
        bond_indices.push_back(b->atom1());
        bond_indices.push_back(b->atom2());
    }
    EXPECT_THAT(bond_indices, UnorderedElementsAre(1, 1, 2, 3));

    EXPECT_THAT(bond_graph.bond(1, 5), IsNull());
    auto bond = bond_graph.bond(2, 1);
    EXPECT_THAT(bond, NotNull());
    EXPECT_EQ(bond, bond_graph.bond(1, 2));
    EXPECT_EQ(bond->atom1(), 1);
    EXPECT_EQ(bond->atom2(), 2);
    EXPECT_EQ(bond->order(), 0);
    EXPECT_FALSE(bond->aromatic());
    EXPECT_TRUE(bond->guessed());
    EXPECT_TRUE(bond->guessed_order());
}

TEST(Bonds, Bond) {
    Bond bond(2, 1);

    EXPECT_EQ(bond.order(), 0);
    bond.set_order(1);
    EXPECT_EQ(bond.order(), 1);

    EXPECT_TRUE(bond.guessed());
    bond.set_guessed(false);
    EXPECT_FALSE(bond.guessed());

    EXPECT_TRUE(bond.guessed_order());
    bond.set_guessed_order(false);
    EXPECT_FALSE(bond.guessed_order());

    EXPECT_FALSE(bond.aromatic());
    bond.set_aromatic(true);
    EXPECT_TRUE(bond.aromatic());

    EXPECT_EQ(bond.atom1(), 1);
    EXPECT_EQ(bond.atom2(), 2);
    bond = Bond(1, 2);
    EXPECT_EQ(bond.atom1(), 1);
    EXPECT_EQ(bond.atom2(), 2);
}
