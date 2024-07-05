#include "utils.hpp"

#include <molpp/internal/Sel.hpp>
#include <molpp/internal/MolData.hpp>
#include <molpp/Atom.hpp>
#include <molpp/AtomSel.hpp>
#include <molpp/ResidueSel.hpp>
#include <molpp/ChainSel.hpp>
#include <molpp/SegmentSel.hpp>
#include <molpp/Error.hpp>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <iterator>

using namespace testing;

//! Type-parameterized test fixture for Sel
template<class SelType>
class SelTest : public ::testing::Test
{
public:
    SelTest()
    : data(create_moldata(4, 1, 4, 4, 2))
    , selection{std::to_array<mol::index_t>({1, 2, 3}), data}
    , selection_as_is{std::to_array<mol::index_t>({2, 2, 0, 1}), data, true}
    {
    }

    mol::internal::MolData data;
    SelType selection;
    SelType selection_as_is;
};

TYPED_TEST_SUITE_P(SelTest);

// Construct from indices and data
TYPED_TEST_P(SelTest, construct_from_indices_and_data)
{
    TypeParam new_selection(std::to_array<mol::index_t>({1, 2}), this->data);

    EXPECT_EQ(new_selection.size(), 2);
    EXPECT_EQ(new_selection.frame(), 0);
    EXPECT_EQ(new_selection[0].index(), 1);
    EXPECT_EQ(new_selection[1].index(), 2);
}

// Construct from non-unique and unsorted indices
TYPED_TEST_P(SelTest, construct_from_non_processed_indices)
{
    TypeParam new_selection(std::to_array<mol::index_t>({2, 1, 1}), this->data, true);

    EXPECT_EQ(new_selection.size(), 3);
    EXPECT_EQ(new_selection.frame(), 0);
    EXPECT_EQ(new_selection[0].index(), 2);
    EXPECT_EQ(new_selection[1].index(), 1);
    EXPECT_EQ(new_selection[2].index(), 1);
}

// Construct from SelIndices and data
TYPED_TEST_P(SelTest, construct_from_selindices_and_data)
{
    mol::internal::SelIndices sel_indices{std::to_array<mol::index_t>({1, 2}), false};
    TypeParam new_selection(sel_indices, this->data);

    EXPECT_EQ(new_selection.size(), 2);
    EXPECT_EQ(new_selection.frame(), 0);
    EXPECT_EQ(new_selection[0].index(), 1);
    EXPECT_EQ(new_selection[1].index(), 2);
}

// Construct only from data
TYPED_TEST_P(SelTest, construct_from_data)
{
    TypeParam new_selection(this->data);

    EXPECT_EQ(new_selection.size(), 4);
    EXPECT_EQ(new_selection.frame(), 0);
    EXPECT_EQ(new_selection[0].index(), 0);
    EXPECT_EQ(new_selection[1].index(), 1);
    EXPECT_EQ(new_selection[2].index(), 2);
    EXPECT_EQ(new_selection[3].index(), 3);
}

// Construct from Entity
TYPED_TEST_P(SelTest, construct_from_entity)
{
    using Entity = mol::internal::SelTraits<TypeParam>::entity_type;
    Entity entity(1, 1, this->data);

    TypeParam new_selection(entity);

    EXPECT_EQ(new_selection.size(), 1);
    EXPECT_EQ(new_selection.frame(), 1);
    EXPECT_EQ(new_selection[0].index(), 1);
}

// Category should be the same as the entity category
TYPED_TEST_P(SelTest, category)
{
    EXPECT_EQ(TypeParam::category(), mol::internal::SelTraits<TypeParam>::entity_type::category());
}

// Frame should be set to the first one if the trajectory has frames
TYPED_TEST_P(SelTest, default_frame_with_trajectory)
{
    EXPECT_EQ(this->selection.frame(), 0);
}

// Frame should be set to nullopt if the trajectory has no frames
TYPED_TEST_P(SelTest, default_frame_without_trajectory)
{
    mol::internal::MolData data_no_frames(1);
    TypeParam selection_no_frames(data_no_frames);

    EXPECT_FALSE(selection_no_frames.frame());
}

// Setting frame to an existing frame should change the frame
TYPED_TEST_P(SelTest, set_valid_frame)
{
    this->selection.set_frame(1);

    EXPECT_EQ(this->selection.frame(), 1);
}

// Setting frame to an invalid frame should throw
TYPED_TEST_P(SelTest, set_invalid_frame)
{
    EXPECT_THROW(this->selection.set_frame(2), mol::Error);
}

// Setting frame to nullopt should update the frame
TYPED_TEST_P(SelTest, set_null_frame)
{
    this->selection.set_frame(std::nullopt);

    EXPECT_FALSE(this->selection.frame());
}

// Size should return the number of selected atoms
TYPED_TEST_P(SelTest, size)
{
    EXPECT_EQ(this->selection.size(), 3);

    EXPECT_EQ(this->selection_as_is.size(), 4);
}

// Sel::contains should return true if the atom is in the selection
TYPED_TEST_P(SelTest, contains)
{
    EXPECT_TRUE(this->selection.contains(1));
    EXPECT_TRUE(this->selection.contains(2));
    EXPECT_TRUE(this->selection.contains(3));

    EXPECT_TRUE(this->selection_as_is.contains(0));
    EXPECT_TRUE(this->selection_as_is.contains(1));
    EXPECT_TRUE(this->selection_as_is.contains(2));
}

// Sel::contains should return false if the atom is not in the selection
TYPED_TEST_P(SelTest, does_not_contain)
{
    EXPECT_FALSE(this->selection.contains(0));
    EXPECT_FALSE(this->selection.contains(4));

    EXPECT_FALSE(this->selection_as_is.contains(3));
}

// Sel::indices should return the selected atom indices
TYPED_TEST_P(SelTest, indices)
{
    EXPECT_THAT(this->selection.indices(), ElementsAre(1, 2, 3));

    EXPECT_THAT(this->selection_as_is.indices(), ElementsAre(2, 2, 0, 1));
}

// Sel::begin should return an iterator to the beginning of the selection
TYPED_TEST_P(SelTest, begin)
{
    auto iter = this->selection.begin();

    EXPECT_EQ((*iter).index(), 1);
}

// Sel::end should return an iterator past to the end of the selection
TYPED_TEST_P(SelTest, end)
{
    auto iter = this->selection.end();

    EXPECT_EQ((*std::prev(iter)).index(), 3);
}

// Sel::begin and Sel::end should form a valid range
TYPED_TEST_P(SelTest, valid_range)
{
    std::vector<mol::index_t> indices;
    for (auto const& entity : this->selection)
    {
        indices.push_back(entity.index());
    }

    EXPECT_THAT(indices, ElementsAre(1, 2, 3));
}

// Sel should be indexable
TYPED_TEST_P(SelTest, indexable)
{
    EXPECT_EQ(this->selection[0].index(), 1);
    EXPECT_EQ(this->selection[1].index(), 2);
    EXPECT_EQ(this->selection[2].index(), 3);

    EXPECT_EQ(this->selection_as_is[0].index(), 2);
    EXPECT_EQ(this->selection_as_is[1].index(), 2);
    EXPECT_EQ(this->selection_as_is[2].index(), 0);
    EXPECT_EQ(this->selection_as_is[3].index(), 1);
}

// Sel::at should return the atom at the given index
TYPED_TEST_P(SelTest, at)
{
    EXPECT_EQ(this->selection.at(0).index(), 1);
    EXPECT_EQ(this->selection.at(1).index(), 2);
    EXPECT_EQ(this->selection.at(2).index(), 3);

    EXPECT_EQ(this->selection_as_is.at(0).index(), 2);
    EXPECT_EQ(this->selection_as_is.at(1).index(), 2);
    EXPECT_EQ(this->selection_as_is.at(2).index(), 0);
    EXPECT_EQ(this->selection_as_is.at(3).index(), 1);
}

// Sel::at should throw if the index is out of bounds
TYPED_TEST_P(SelTest, at_out_of_bounds)
{
    EXPECT_THROW(this->selection.at(3), mol::Error);

    EXPECT_THROW(this->selection_as_is.at(4), mol::Error);
}

// Sel::by_index should return the atom whose index is the given one
TYPED_TEST_P(SelTest, by_index)
{
    EXPECT_EQ(this->selection.by_index(1).index(), 1);
    EXPECT_EQ(this->selection.by_index(2).index(), 2);
    EXPECT_EQ(this->selection.by_index(3).index(), 3);

    EXPECT_EQ(this->selection_as_is.by_index(0).index(), 0);
    EXPECT_EQ(this->selection_as_is.by_index(1).index(), 1);
    EXPECT_EQ(this->selection_as_is.by_index(2).index(), 2);
}

// Sel::by_index should throw if the index is not in the selection
TYPED_TEST_P(SelTest, by_index_not_found)
{
    EXPECT_THROW(this->selection.by_index(0), mol::Error);
    EXPECT_THROW(this->selection.by_index(4), mol::Error);

    EXPECT_THROW(this->selection_as_is.by_index(3), mol::Error);
}

// Iterator should be dereferenceable
TYPED_TEST_P(SelTest, iterator_dereferenceable)
{
    auto iter = this->selection.begin();

    EXPECT_EQ((*iter).index(), 1);
}

// Iterator should be pre-incrementable
TYPED_TEST_P(SelTest, iterator_pre_incrementable)
{
    auto iter = this->selection.begin();
    ++iter;

    EXPECT_EQ((*iter).index(), 2);
}

// Iterator should be post-incrementable
TYPED_TEST_P(SelTest, iterator_post_incrementable)
{
    auto iter = this->selection.begin();
    auto previous = iter++;

    EXPECT_EQ((*previous).index(), 1);
    EXPECT_EQ((*iter).index(), 2);
}

// Iterator should be pre-decrementable
TYPED_TEST_P(SelTest, iterator_pre_decrementable)
{
    auto iter = this->selection.end();
    --iter;

    EXPECT_EQ((*iter).index(), 3);
}

// Iterator should be post-decrementable
TYPED_TEST_P(SelTest, iterator_post_decrementable)
{
    auto iter = this->selection.end();
    auto previous = iter--;

    EXPECT_EQ(previous, this->selection.end());
    EXPECT_EQ((*iter).index(), 3);
}

// Iterator should be incrementable by an integer
TYPED_TEST_P(SelTest, iterator_incrementable_by_integer)
{
    auto iter = this->selection.begin();
    iter += 1;

    EXPECT_EQ((*iter).index(), 2);
}

// Iterator should be decrementable by an integer
TYPED_TEST_P(SelTest, iterator_decrementable_by_integer)
{
    auto iter = this->selection.end();
    iter -= 1;

    EXPECT_EQ((*iter).index(), 3);
}

// Iterator should be subtractable
TYPED_TEST_P(SelTest, iterator_subtractable)
{
    auto iter1 = this->selection.begin();
    auto iter2 = this->selection.end();

    EXPECT_EQ(iter2 - iter1, 3);
}

// Iterator should be comparable
TYPED_TEST_P(SelTest, iterator_comparable)
{
    auto iter1 = this->selection.begin();
    auto iter2 = this->selection.begin();

    EXPECT_EQ(iter1, iter2);
}

// Iterator should be inequality-comparable
TYPED_TEST_P(SelTest, iterator_inequality_comparable)
{
    auto iter1 = this->selection.begin();
    auto iter2 = this->selection.end();

    EXPECT_NE(iter1, iter2);
}

REGISTER_TYPED_TEST_SUITE_P(SelTest, construct_from_indices_and_data, construct_from_non_processed_indices, construct_from_selindices_and_data, construct_from_data, construct_from_entity, category, default_frame_with_trajectory, default_frame_without_trajectory, set_valid_frame, set_invalid_frame, set_null_frame, size, contains, does_not_contain, indices, begin, end, valid_range, indexable, at, at_out_of_bounds, by_index, by_index_not_found, iterator_dereferenceable, iterator_pre_incrementable, iterator_post_incrementable, iterator_pre_decrementable, iterator_post_decrementable, iterator_incrementable_by_integer, iterator_decrementable_by_integer, iterator_subtractable, iterator_comparable, iterator_inequality_comparable);

class SelMock;

namespace mol::internal
{

template<>
struct SelTraits<SelMock>
{
    using entity_type = Atom;
};

} // namespace mol::internal

//! Class that implements the Sel interface for testing
struct SelMock : public mol::internal::Sel<SelMock>
{
    using entity_type = mol::Atom;
    using mol::internal::Sel<SelMock>::Sel;

    template<class, class>
    friend class mol::internal::Sel;
};

// Test all types that implement the Sel interface
using SelTypes = ::testing::Types<SelMock, mol::AtomSel, mol::ResidueSel, mol::ChainSel, mol::SegmentSel>;
INSTANTIATE_TYPED_TEST_SUITE_P(SelInterface, SelTest, SelTypes);
