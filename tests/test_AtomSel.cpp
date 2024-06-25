#include "auxiliary.hpp"

#include <molpp/internal/MolData.hpp>
#include <molpp/AtomSel.hpp>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <array>
#include <utility>

using namespace testing;

//! Test fixture for AtomSel
class AtomSelTest : public ::testing::Test
{
public:
    AtomSelTest()
    : data(create_moldata(4, 1, 4, 4, 2))
    , selection{std::to_array<mol::index_t>({1, 2, 3}), data}
    , const_selection{std::to_array<mol::index_t>({1, 2, 3}), data}
    {
    }

    mol::internal::MolData data;
    mol::AtomSel selection;
    mol::AtomSel const const_selection;
};

// AtomSel::positions should the positions of only the selected atoms
TEST_F(AtomSelTest, positions)
{
     mol::Coord3 positions = selection.positions();
     mol::Coord3 const const_positions = const_selection.positions();

    EXPECT_THAT(positions.reshaped(), ElementsAre(1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0));
    EXPECT_THAT(const_positions.reshaped(), ElementsAre(1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0));
}

// AtomSel::bonds should return a selection with the atoms bonded to the selected atoms
TEST_F(AtomSelTest, bonded)
{
    mol::AtomSel bonded = selection.bonded();

    EXPECT_THAT(bonded.indices(), ElementsAre(0, 1, 2, 3));
}

// AtomSel::bonds should return the bonds to the selected atoms
TEST_F(AtomSelTest, bonds)
{
    std::vector<std::pair<mol::index_t, mol::index_t>> bonds_indices;
    for (std::shared_ptr<mol::Bond> bond : selection.bonds())
    {
        bonds_indices.push_back(std::make_pair<mol::index_t, mol::index_t>(bond->atom1(), bond->atom2()));
    }

    EXPECT_THAT(bonds_indices, UnorderedElementsAre(Pair(0, 1), Pair(1, 2), Pair(2, 3)));
}
