#include "utils.hpp"

#include <molpp/internal/MolData.hpp>
#include <molpp/Segment.hpp>
#include <molpp/Residue.hpp>
#include <molpp/ResidueSel.hpp>
#include <molpp/Error.hpp>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace mol;
using namespace testing;

//! Test fixture for Segment class
class SegmentTest : public ::testing::Test
{
public:
    SegmentTest()
    : data(create_moldata(3, 1, 1, 3, 1))
    , segment(1, 0, data)
    , const_segment(1, 0, data)
    , null_frame_segment(0, std::nullopt, data)
    {}

    mol::internal::MolData data;
    Segment segment;
    Segment const const_segment;
    Segment null_frame_segment;
};

// Segment::category should return the correct category
TEST_F(SegmentTest, category)
{
    EXPECT_EQ(Segment::category(), EntityCategory::Segment);
}

// Segment::size should return the number of residues in the segment
TEST_F(SegmentTest, size)
{
    EXPECT_EQ(segment.size(), 1);
}

// Segment::index should return the correct index
TEST_F(SegmentTest, index)
{
    EXPECT_EQ(segment.index(), 1);
    EXPECT_EQ(const_segment.index(), 1);
}

// Segment::indices should return a list with only the index
TEST_F(SegmentTest, indices)
{
    EXPECT_THAT(segment.indices(), ElementsAre(1));
    EXPECT_THAT(const_segment.indices(), ElementsAre(1));
}

// Segment::frame should return the correct frame
TEST_F(SegmentTest, frames)
{
    EXPECT_EQ(segment.frame(), 0);
    EXPECT_EQ(const_segment.frame(), 0);
    EXPECT_FALSE(null_frame_segment.frame());
}

// Setting a valid frame should update the frame
TEST_F(SegmentTest, set_valid_frame)
{
    segment.set_frame(0);
    EXPECT_EQ(segment.frame(), 0);
}

// Setting a frame to nullopt should update the frame
TEST_F(SegmentTest, set_null_frame)
{
    segment.set_frame(std::nullopt);
    EXPECT_FALSE(segment.frame());
}

// Setting an invalid frame should throw a Error
TEST_F(SegmentTest, set_invalid_frame)
{
    EXPECT_THROW(segment.set_frame(1), Error);
}

// Segment::name should return the correct name
TEST_F(SegmentTest, name)
{
    EXPECT_EQ(segment.name(), "B");
    EXPECT_EQ(const_segment.name(), "B");
}

// Segment::set_name should update the name
TEST_F(SegmentTest, set_name)
{
    segment.set_name("A");

    EXPECT_EQ(segment.name(), "A");
}

// Segment::add_residue with an index should make the residue exclusive to the new segment
TEST_F(SegmentTest, add_residue_from_index)
{
    Residue new_residue(0, std::nullopt, data);
    Segment old_segment = new_residue.segment().value();
    ASSERT_NE(old_segment, segment);

    segment.add_residue(new_residue.index());
    ResidueSel segment_residues(segment);
    ResidueSel old_segment_residues(old_segment);

    EXPECT_THAT(segment_residues.indices(), UnorderedElementsAre(0, 1));
    EXPECT_THAT(old_segment_residues.indices(), ElementsAre());
    EXPECT_EQ(new_residue.segment_index(), segment.index());
    EXPECT_EQ(segment.size(), 2);
    EXPECT_EQ(old_segment.size(), 0);
}

// Segment::add_residue with an Residue should make the residue exclusive to the new segment
TEST_F(SegmentTest, add_residue_from_residue)
{
    Residue new_residue(0, std::nullopt, data);
    Segment old_segment = new_residue.segment().value();
    ASSERT_NE(old_segment, segment);

    segment.add_residue(new_residue);
    ResidueSel segment_residues(segment);
    ResidueSel old_segment_residues(old_segment);

    EXPECT_THAT(segment_residues.indices(), UnorderedElementsAre(0, 1));
    EXPECT_THAT(old_segment_residues.indices(), ElementsAre());
    EXPECT_EQ(new_residue.segment_index(), segment.index());
    EXPECT_EQ(segment.size(), 2);
    EXPECT_EQ(old_segment.size(), 0);
}
