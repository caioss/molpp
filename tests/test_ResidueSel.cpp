#include "auxiliary.hpp"

#include <molpp/internal/MolData.hpp>
#include <molpp/ResidueSel.hpp>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <array>
#include <utility>

using namespace testing;

//! Test fixture for ResidueSel
class ResidueSelTest : public ::testing::Test
{
public:
    ResidueSelTest()
    : data(create_moldata(4, 1, 4, 4, 2))
    , selection{std::to_array<mol::index_t>({1, 2, 3}), data}
    , const_selection{std::to_array<mol::index_t>({1, 2, 3}), data}
    {
    }

    mol::internal::MolData data;
    mol::ResidueSel selection;
    mol::ResidueSel const const_selection;
};

// ResidueSel::as_atom_indices should return the indices of the selected atoms
TEST_F(ResidueSelTest, as_atom_indices)
{
    EXPECT_THAT(selection.as_atom_indices(), ElementsAre(1, 2, 3));
}

// ResidueSel::from_atom_indices should return the atom indices
TEST_F(ResidueSelTest, from_atom_indices)
{
    mol::ResidueSel::indices_type const indices{1, 2};
    mol::ResidueSel::indices_type from_atom_indices = mol::ResidueSel::from_atom_indices(indices, data);

    EXPECT_THAT(from_atom_indices, ElementsAre(1, 2));
}
