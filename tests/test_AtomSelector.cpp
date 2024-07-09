#include "utils.hpp"

#include <molpp/AtomSelector.hpp>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace mol;
using namespace mol::internal;
using namespace testing;

class AtomSelectorTest : public testing::Test
{
public:
    AtomSelectorTest()
    : data(create_moldata(4, 1, 4, 4, 2))
    , selector("resid 1", data)
    {}

    mol::internal::MolData data;
    AtomSelector selector;
};

// Construct an AtomSelector with a valid selection
TEST_F(AtomSelectorTest, construct_valid_selection)
{
    AtomSelector new_selector("resid 1", data);
}

// Construct an AtomSelector with an invalid selection
TEST_F(AtomSelectorTest, construct_invalid_selection)
{
    EXPECT_THROW(AtomSelector("this is absolutely invalid", data), mol::Error);
}

// Apply a non-empty selection should return a non-empty AtomSel regarding the frame
TEST_F(AtomSelectorTest, apply_non_empty_selection)
{
    AtomSel atoms_frame_0 = selector.apply(0);
    AtomSel atoms_frame_1 = selector.apply(1);
    AtomSel atoms_null_frame = selector.apply(std::nullopt);

    EXPECT_THAT(atoms_frame_0.indices(), ElementsAre(1));
    EXPECT_THAT(atoms_frame_1.indices(), ElementsAre(1));
    EXPECT_THAT(atoms_null_frame.indices(), ElementsAre(1));
}

// Apply an empty selection should return an empty AtomSel
TEST_F(AtomSelectorTest, apply_empty_selection)
{
    AtomSelector empty_selector("resid 999", data);
    AtomSel atoms = empty_selector.apply(0);

    EXPECT_THAT(atoms.indices(), IsEmpty());
}

// Apply should return a selection with the correct frame
TEST_F(AtomSelectorTest, apply_frame)
{
    AtomSel atoms_frame_0 = selector.apply(0);
    AtomSel atoms_frame_1 = selector.apply(1);
    AtomSel atoms_null_frame = selector.apply(std::nullopt);

    EXPECT_EQ(atoms_frame_0.frame(), 0);
    EXPECT_EQ(atoms_frame_1.frame(), 1);
    EXPECT_FALSE(atoms_null_frame.frame());
}

// Apply with an invalid frame should throw
TEST_F(AtomSelectorTest, apply_invalid_frame)
{
    EXPECT_THROW(selector.apply(2), mol::Error);
}

struct SelectionTestData
{
    std::string selection;
    std::vector<index_t> expected;
};

class SelectionTest : public testing::TestWithParam<SelectionTestData>
{
public:
    SelectionTest()
    : data(create_moldata(4, 1, 4, 4, 2))
    , selector(GetParam().selection, data)
    {}

    mol::internal::MolData data;
    AtomSelector selector;
};

// AtomSelector should select correct atoms for a given selection string
TEST_P(SelectionTest, select_atoms)
{
    AtomSel atoms = selector.apply(0);

    EXPECT_THAT(atoms.indices(), ElementsAreArray(GetParam().expected));
}

// Resid property
INSTANTIATE_TEST_SUITE_P(Resid, SelectionTest, Values(
    SelectionTestData{"resid 0", {0}},
    SelectionTestData{"resid 1", {1}},
    SelectionTestData{"resid 2", {2}},
    SelectionTestData{"resid 3", {3}},
    SelectionTestData{"resid 4", {}},
    SelectionTestData{"resid 1:3", {1, 2, 3}},
    SelectionTestData{"resid 0 1:3", {0, 1, 2, 3}},
    SelectionTestData{"resid 0:2 3", {0, 1, 2, 3}},
    SelectionTestData{"resid 0 2 3", {0, 2, 3}},
    SelectionTestData{"resid 0 1 2:3", {0, 1, 2, 3}},
    SelectionTestData{"resid 0:1 2 3", {0, 1, 2, 3}}
));

// Boolean operators
INSTANTIATE_TEST_SUITE_P(Boolean, SelectionTest, Values(
    SelectionTestData{"resid 1 or resid 2", {1, 2}},
    SelectionTestData{"resid 1 and resid 2", {}},
    SelectionTestData{"resid 1 or resid 2 or resid 3", {1, 2, 3}},
    SelectionTestData{"resid 1 and resid 2 or resid 3", {3}},
    SelectionTestData{"resid 1 or resid 2 and resid 3", {}},
    SelectionTestData{"resid 1 and resid 2 and resid 3", {}},
    SelectionTestData{"resid 0 or resid 1 or resid 2 or resid 3", {0, 1, 2, 3}},
    SelectionTestData{"resid 0 and resid 1 or resid 2 or resid 3", {2, 3}},
    SelectionTestData{"resid 0 or resid 1 and resid 2 or resid 3", {3}},
    SelectionTestData{"resid 0 and resid 1 and resid 2 or resid 3", {3}},
    SelectionTestData{"resid 0 or resid 1 or resid 2 and resid 3", {}},
    SelectionTestData{"resid 0 and resid 1 or resid 2 and resid 3", {}},
    SelectionTestData{"resid 0 or resid 1 and resid 2 and resid 3", {}},
    SelectionTestData{"resid 0 and resid 1 and resid 2 and resid 3", {}},
    SelectionTestData{"not resid 0", {1, 2, 3}},
    SelectionTestData{"not resid 0 or resid 1", {1, 2, 3}},
    SelectionTestData{"not resid 0 or resid 0", {0, 1, 2, 3}},
    SelectionTestData{"not resid 0 and resid 1", {1}},
    SelectionTestData{"not resid 0 and resid 0", {}},
    SelectionTestData{"resid 0 or not resid 1", {0, 2, 3}},
    SelectionTestData{"resid 0 or not resid 0", {0, 1, 2, 3}},
    SelectionTestData{"resid 0 and not resid 1", {0}},
    SelectionTestData{"resid 0 and not resid 0", {}}
));

// Parentheses
INSTANTIATE_TEST_SUITE_P(Parentheses, SelectionTest, Values(
    SelectionTestData{"(resid 0)", {0}},
    SelectionTestData{"(resid 0 or resid 1)", {0, 1}},
    SelectionTestData{"(resid 0 or resid 1) or resid 2", {0, 1, 2}},
    SelectionTestData{"(resid 0 or resid 1) and resid 2", {}},
    SelectionTestData{"(resid 0 or resid 1) or (resid 2 or resid 3)", {0, 1, 2, 3}},
    SelectionTestData{"(resid 0 or resid 1) and (resid 2 or resid 3)", {}},
    SelectionTestData{"(resid 0 and resid 1) or (resid 2 and resid 3)", {}},
    SelectionTestData{"(resid 0 and resid 1) and (resid 2 or resid 3)", {}},
    SelectionTestData{"not (resid 0)", {1, 2, 3}},
    SelectionTestData{"not (resid 0 or resid 1)", {2, 3}},
    SelectionTestData{"(resid 0 or resid 1) or not resid 2", {0, 1, 3}},
    SelectionTestData{"not (resid 0 or resid 1) and resid 2", {2}},
    SelectionTestData{"(resid 0 or resid 1) or not (resid 2 or resid 3)", {0, 1}},
    SelectionTestData{"(not resid 0 or resid 1) and not (resid 2 or resid 3)", {1}},
    SelectionTestData{"(resid 0 and not resid 1) or (resid 2 and resid 3)", {0}},
    SelectionTestData{"(resid 0 and resid 1) and (not resid 2 or resid 3)", {}},
    SelectionTestData{"not ((resid 0 and resid 1) and (resid 2 or resid 3))", {0, 1, 2, 3}}
));

// All selector
INSTANTIATE_TEST_SUITE_P(All, SelectionTest, Values(
    SelectionTestData{"all", {0, 1, 2, 3}},
    SelectionTestData{"not all", {}},
    SelectionTestData{"resid 1 and all", {1}},
    SelectionTestData{"resid 1 and not all", {}},
    SelectionTestData{"resid 1 or all", {0, 1, 2, 3}}
));
