#include "utils.hpp"

#include "selections/SelectionStack.hpp"
#include "selections/properties.hpp"
#include <molpp/internal/MolData.hpp>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace mol;
using namespace mol::internal;
using namespace testing;

//! Test fixture for ResidSelection selection nodes.
class ResidSelectionTest : public testing::Test
{
public:
    ResidSelectionTest()
    : data(create_moldata(5, 1, 1, 1, 1))
    {
        indices.available->insert(0);
        indices.available->insert(1);
        indices.available->insert(2);
        indices.available->insert(3);

        stack.push_indices(indices);

        resids.add_number({"0"});
        resids.add_range({"2", "4"});
    }

    MolData data;
    SelectionIndices indices;
    SelectionStack stack;
    NumberSet resids;
};

// ResidSelection with an empty residue set should not select any atoms.
TEST_F(ResidSelectionTest, empty_resid_set)
{
    ResidSelection selection(NumberSet{});

    selection.evaluate(stack, data, 0);

    EXPECT_THAT(*indices.selected, IsEmpty());
}

// ResidSelection should select atoms that are in the given residue set and available to be selected.
TEST_F(ResidSelectionTest, from_resid_set)
{
    ResidSelection selection(std::move(resids));

    selection.evaluate(stack, data, 0);

    EXPECT_THAT(*indices.selected, UnorderedElementsAre(0, 2, 3));
}

// ResidSelection::evaluate_atom should return true if the atom is in the residue set, and false otherwise.
TEST_F(ResidSelectionTest, evaluate_atom)
{
    ResidSelection selection(std::move(resids));

    EXPECT_TRUE(selection.evaluate_atom(0, data));
    EXPECT_FALSE(selection.evaluate_atom(1, data));
    EXPECT_TRUE(selection.evaluate_atom(2, data));
    EXPECT_TRUE(selection.evaluate_atom(3, data));
    EXPECT_TRUE(selection.evaluate_atom(4, data));
}
