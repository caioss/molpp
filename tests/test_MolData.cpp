#include <molpp/Atom.hpp>
#include <molpp/Residue.hpp>
#include <molpp/internal/MolData.hpp>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace testing;
using namespace mol;
using namespace mol::internal;

class MolDataTest : public ::testing::Test
{
protected:
    MolDataTest()
    : size{3}
    , num_frames{5}
    , data{size}
    , const_data{data}
    {
        data.add_entity<Atom>(size);
        data.add_property<Atom, Name>(true);
        data.add_property<Atom, Charge>(false);

        for (size_t i = 0; i < num_frames; i++)
        {
            data.add_frame();
        }
    }

    size_t const size;
    size_t const num_frames;
    MolData data;
    MolData const& const_data;
};

TEST_F(MolDataTest, NotRegisteredEntitySize)
{
    EXPECT_THROW(data.entity_size<Residue>(), mol::MolError);
}

TEST_F(MolDataTest, AddEntity)
{
    bool const result = data.add_entity<Residue>(size);

    EXPECT_TRUE(result);
    EXPECT_EQ(const_data.entity_size<Residue>(), size);
}

TEST_F(MolDataTest, AddEntityTwice)
{
    data.add_entity<Residue>(size);
    bool const result = data.add_entity<Residue>(size);

    EXPECT_FALSE(result);
}

TEST_F(MolDataTest, ResizeEntity)
{
    size_t const new_size = size + 2;
    data.set_entity_size<Atom>(new_size);

    EXPECT_EQ(const_data.entity_size<Atom>(), new_size);
    for (size_t frame = 0; frame < num_frames; frame++)
    {
        Name* name = data.property_at<Atom, Name>(frame);
        ASSERT_THAT(name, NotNull()) << "Frame " << frame;
        EXPECT_EQ(name->size(), new_size) << "Frame " << frame;
    }
}

TEST_F(MolDataTest, ResizeNotRegisteredEntity)
{
    EXPECT_THROW(data.set_entity_size<Residue>(3), mol::MolError);
}

TEST_F(MolDataTest, DefaultNumberOfFrames)
{
    MolData empty_data(0);
    EXPECT_EQ(empty_data.num_frames(), 0);
}

TEST_F(MolDataTest, AddFrames)
{
    for (size_t i = num_frames; i < num_frames + 3; i++)
    {
        Frame frame = data.add_frame();
        Name* name = data.property_at<Atom, Name>(frame);

        EXPECT_EQ(frame, i) << "Frame " << i;
        EXPECT_EQ(data.num_frames(), i + 1) << "Frame " << i;
        EXPECT_THAT(name, NotNull()) << "Frame " << i;
    }
}

TEST_F(MolDataTest, RemoveFrames)
{
    data.remove_frame(1);
    data.remove_frame(2);

    EXPECT_EQ(data.num_frames(), num_frames - 2);
    for (size_t frame = 0; frame < num_frames - 2; frame++)
    {
        Name* name = data.property_at<Atom, Name>(frame);
        EXPECT_THAT(name, NotNull()) << "Frame " << frame;
    }
    Name* name = data.property_at<Atom, Name>(num_frames - 2);
    EXPECT_THAT(name, IsNull());
}

TEST_F(MolDataTest, GetTimeBasedProperty)
{
    for (size_t frame = 0; frame < num_frames; frame++)
    {
        Name* name = data.property_at<Atom, Name>(frame);
        Name const* const_name = const_data.property_at<Atom, Name>(frame);
        EXPECT_THAT(name, NotNull()) << "Frame " << frame;
        EXPECT_THAT(const_name, NotNull()) << "Frame " << frame;
        EXPECT_EQ(name, const_name) << "Frame " << frame;
    }

    Name* name = data.property_at<Atom, Name>(num_frames);
    Name const* const_name = const_data.property_at<Atom, Name>(num_frames);
    EXPECT_THAT(name, IsNull());
    EXPECT_THAT(const_name, IsNull());
}

TEST_F(MolDataTest, GetNonTimeBasedProperty)
{
    Charge* charge = data.property_at<Atom, Charge>(0);

    for (size_t frame = 0; frame < num_frames; frame++)
    {
        Charge* frame_charge = data.property_at<Atom, Charge>(frame);

        EXPECT_EQ(charge, frame_charge) << "Frame " << frame;
    }
}

TEST_F(MolDataTest, AddTimeBasedProperty)
{
    ResName* resname = data.add_property<Atom, ResName>(true);

    EXPECT_THAT(resname, NotNull());

    for (size_t frame = 0; frame < num_frames; frame++)
    {
        ResName* frame_resname = data.property_at<Atom, ResName>(frame);

        EXPECT_THAT(frame_resname, NotNull()) << "Frame " << frame;
    }
}

TEST_F(MolDataTest, AddNonTimeBasedProperty)
{
    ResName* resname = data.add_property<Atom, ResName>(false);

    EXPECT_THAT(resname, NotNull());

    for (size_t frame = 0; frame < num_frames; frame++)
    {
        ResName* frame_resname = data.property_at<Atom, ResName>(frame);

        EXPECT_EQ(frame_resname, resname) << "Frame " << frame;
    }
}

TEST_F(MolDataTest, AddPropertyOnNotRegisteredEntity)
{
    EXPECT_THROW((data.add_property<Residue, ResName>(true)), mol::MolError);
}
