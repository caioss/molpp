#include "selections/SelectionStack.hpp"
#include "selections/PropertySelection.hpp"
#include <molpp/internal/MolData.hpp>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace mol;
using namespace mol::internal;
using namespace testing;

class HelperSelection : public PropertySelection<HelperSelection>
{
public:
    HelperSelection(bool accept)
    : m_accept{accept}
    {}

    bool evaluate_atom(index_t const /*atom_index*/, MolData const& /*data*/) const
    {
        return m_accept;
    }

private:
    bool m_accept;
};

//! Test fixture for PropertySelection selection nodes.
class PropertySelectionTest : public testing::Test
{
public:
    PropertySelectionTest()
    : data(0)
    {
        indices.available->insert(0);
        indices.available->insert(1);
        indices.available->insert(2);

        stack.push_indices(indices);
    }

    MolData data;
    SelectionIndices indices;
    SelectionStack stack;
};

// Atom index should be selected if evaluate_atom returns true.
TEST_F(PropertySelectionTest, select_atom)
{
    HelperSelection selection(true);
    selection.evaluate(stack, data, 0);

    EXPECT_THAT(*indices.selected, UnorderedElementsAre(0, 1, 2));
}

// Atom index should not be selected if evaluate_atom returns false.
TEST_F(PropertySelectionTest, deselect_atom)
{
    HelperSelection selection(false);
    selection.evaluate(stack, data, 0);

    EXPECT_TRUE(indices.selected->empty());
}
