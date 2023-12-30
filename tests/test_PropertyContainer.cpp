#include <molpp/internal/PropertyContainer.hpp>
#include <molpp/Atom.hpp>
#include <molpp/Residue.hpp>
#include <molpp/Property.hpp>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace testing;
using namespace mol;
using namespace mol::internal;

class PropertyContainerTest : public ::testing::Test
{
protected:
    PropertyContainerTest()
    {}

    PropertyContainer container;
};

TEST_F(PropertyContainerTest, DefaultSize)
{
    container.add<Atom, Mass>(false);
    container.add<Residue, Charge>(false);
    Mass* mass = container.get<Atom, Mass>(0);
    Charge* charge = container.get<Residue, Charge>(0);

    EXPECT_EQ(container.size<Atom>(), 0);
    EXPECT_EQ(container.size<Residue>(), 0);
    EXPECT_EQ(mass->size(), 0);
    EXPECT_EQ(charge->size(), 0);
}

TEST_F(PropertyContainerTest, SetSizeOnMultipleProperties)
{
    Mass* atom_mass = container.add<Atom, Mass>(false);
    Charge* atom_charge = container.add<Atom, Charge>(false);
    Mass* residue_mass = container.add<Residue, Mass>(false);
    Charge* residue_charge = container.add<Residue, Charge>(false);

    container.set_size<Atom>(2);
    container.set_size<Residue>(3);

    EXPECT_EQ(container.size<Atom>(), 2);
    EXPECT_EQ(container.size<Residue>(), 3);

    EXPECT_EQ(atom_mass->size(), 2);
    EXPECT_EQ(atom_charge->size(), 2);
    EXPECT_EQ(residue_mass->size(), 3);
    EXPECT_EQ(residue_charge->size(), 3);
}

TEST_F(PropertyContainerTest, SetSizeBeforeAddFrame)
{
    container.add<Atom, Mass>(true);
    container.add<Residue, Charge>(true);

    container.set_size<Atom>(2);
    container.set_size<Residue>(3);

    size_t const num_frames = 3;

    for (size_t frame = 0; frame < num_frames; frame++)
    {
        container.add_frame();
    }

    EXPECT_EQ(container.size<Atom>(), 2);
    EXPECT_EQ(container.size<Residue>(), 3);

    for (size_t frame = 0; frame < num_frames; frame++)
    {
        Mass* mass = container.get<Atom, Mass>(frame);
        Charge* charge = container.get<Residue, Charge>(frame);

        EXPECT_EQ(mass->size(), 2) << frame;
        EXPECT_EQ(charge->size(), 3) << frame;
    }
}

TEST_F(PropertyContainerTest, SetSizeAfterAddFrame)
{
    container.add<Atom, Mass>(true);
    container.add<Residue, Charge>(true);

    size_t const num_frames = 3;

    for (size_t frame = 0; frame < num_frames; frame++)
    {
        container.add_frame();
    }

    container.set_size<Atom>(2);
    container.set_size<Residue>(3);

    EXPECT_EQ(container.size<Atom>(), 2);
    EXPECT_EQ(container.size<Residue>(), 3);

    for (size_t frame = 0; frame < num_frames; frame++)
    {
        Mass* mass = container.get<Atom, Mass>(frame);
        Charge* charge = container.get<Residue, Charge>(frame);

        EXPECT_EQ(mass->size(), 2) << frame;
        EXPECT_EQ(charge->size(), 3) << frame;
    }
}

TEST_F(PropertyContainerTest, SetSizeBeforeAddFrameNotTimeBased)
{
    Mass* default_mass = container.add<Atom, Mass>(false);
    Charge* default_charge = container.add<Residue, Charge>(false);

    container.set_size<Atom>(2);
    container.set_size<Residue>(3);

    size_t const num_frames = 3;

    for (size_t frame = 0; frame < num_frames; frame++)
    {
        container.add_frame();
    }

    EXPECT_EQ(container.size<Atom>(), 2);
    EXPECT_EQ(container.size<Residue>(), 3);

    for (size_t frame = 0; frame < num_frames; frame++)
    {
        Mass* mass = container.get<Atom, Mass>(frame);
        Charge* charge = container.get<Residue, Charge>(frame);

        EXPECT_EQ(default_mass, mass) << frame;
        EXPECT_EQ(default_charge, charge) << frame;
        EXPECT_EQ(mass->size(), 2) << frame;
        EXPECT_EQ(charge->size(), 3) << frame;
    }
}

TEST_F(PropertyContainerTest, SetSizeAfterAddFrameNotTimeBased)
{
    Mass* default_mass = container.add<Atom, Mass>(false);
    Charge* default_charge = container.add<Residue, Charge>(false);

    size_t const num_frames = 3;

    for (size_t frame = 0; frame < num_frames; frame++)
    {
        container.add_frame();
    }

    container.set_size<Atom>(2);
    container.set_size<Residue>(3);

    EXPECT_EQ(container.size<Atom>(), 2);
    EXPECT_EQ(container.size<Residue>(), 3);

    for (size_t frame = 0; frame < num_frames; frame++)
    {
        Mass* mass = container.get<Atom, Mass>(frame);
        Charge* charge = container.get<Residue, Charge>(frame);

        EXPECT_EQ(default_mass, mass) << frame;
        EXPECT_EQ(default_charge, charge) << frame;
        EXPECT_EQ(mass->size(), 2) << frame;
        EXPECT_EQ(charge->size(), 3) << frame;
    }
}

TEST_F(PropertyContainerTest, GetNotTimeBasedWithoutFrames)
{
    Mass* from_add = container.add<Atom, Mass>(false);
    ASSERT_THAT(from_add, NotNull());

    for (size_t i = 0; i < 3; i++)
    {
        Mass* from_get = container.get<Atom, Mass>(i);

        EXPECT_THAT(from_get, NotNull()) << i;
        EXPECT_EQ(from_add, from_get) << i;
    }
}

TEST_F(PropertyContainerTest, GetTimeBasedWithoutFrames)
{
    Mass* from_add = container.add<Atom, Mass>(true);
    ASSERT_THAT(from_add, IsNull());

    for (size_t i = 0; i < 3; i++)
    {
        Mass* from_get = container.get<Atom, Mass>(i);

        EXPECT_THAT(from_get, IsNull()) << i;
    }
}

TEST_F(PropertyContainerTest, GetNotTimeBasedWithFrames)
{
    size_t const num_frames = 3;
    for (size_t frame = 0; frame < num_frames; frame++)
    {
        container.add_frame();
    }

    Mass* default_mass = container.add<Atom, Mass>(false);

    for (size_t i = 0; i < num_frames; i++)
    {
        Mass* mass = container.get<Atom, Mass>(i);
        EXPECT_THAT(mass, NotNull()) << i;
        EXPECT_EQ(mass, default_mass) << i;
    }
}

TEST_F(PropertyContainerTest, GetTimeBasedWithFrames)
{
    container.add<Atom, Mass>(true);

    size_t const num_frames = 3;
    for (size_t frame = 0; frame < num_frames; frame++)
    {
        container.add_frame();
    }

    for (size_t i = 0; i < num_frames; i++)
    {
        Mass* mass_i = container.get<Atom, Mass>(i);
        EXPECT_THAT(mass_i, NotNull()) << i;

        for (size_t j = i + 1; j < num_frames; j++)
        {
            Mass* mass_j = container.get<Atom, Mass>(j);

            EXPECT_THAT(mass_j, NotNull()) << j;
            EXPECT_NE(mass_i, mass_j) << i << "/" << j;
        }
    }
}

TEST_F(PropertyContainerTest, GetNotRegistered)
{
    container.add<Atom, Mass>(false);

    for (size_t i = 0; i < 3; i++)
    {
        Charge* charge = container.get<Atom, Charge>(i);

        EXPECT_THAT(charge, IsNull()) << i;
    }
}

TEST_F(PropertyContainerTest, AddBeforeAddFrame)
{
    Mass* time_based = container.add<Atom, Mass>(true);
    Charge* not_time_based = container.add<Atom, Charge>(false);

    container.add_frame();

    EXPECT_THAT(time_based, IsNull());
    EXPECT_THAT(not_time_based, NotNull());
}

TEST_F(PropertyContainerTest, AddAfterAddFrame)
{
    container.add_frame();

    Mass* time_based = container.add<Atom, Mass>(true);
    Charge* not_time_based = container.add<Atom, Charge>(false);

    EXPECT_THAT(time_based, NotNull());
    EXPECT_THAT(not_time_based, NotNull());
}

TEST_F(PropertyContainerTest, AddMultipleTimes)
{
    Mass* first = container.add<Atom, Mass>(false);
    Mass* second = container.add<Atom, Mass>(false);
    Mass* third = container.add<Atom, Mass>(true);

    EXPECT_THAT(first, NotNull());
    EXPECT_THAT(second, NotNull());
    EXPECT_THAT(third, NotNull());
    EXPECT_EQ(first, second);
    EXPECT_EQ(first, third);
}

TEST_F(PropertyContainerTest, DefaultNumberOfFrames)
{
    EXPECT_EQ(container.num_frames(), 0);
}

TEST_F(PropertyContainerTest, AddFrames)
{
    for (size_t frame = 0; frame < 3; frame++)
    {
        EXPECT_EQ(container.add_frame(), frame) << frame;
        EXPECT_EQ(container.num_frames(), frame + 1) << frame;
    }
}

TEST_F(PropertyContainerTest, AddFrameKeepSize)
{
    container.add<Atom, Mass>(true);
    container.add<Atom, Charge>(true);
    container.set_size<Atom>(3);

    size_t const num_frames = 3;
    for (size_t frame = 0; frame < num_frames; frame++)
    {
        container.add_frame();
    }

    for (size_t frame = 0; frame < num_frames; frame++)
    {
        Mass* mass = container.get<Atom, Mass>(frame);
        Charge* charge = container.get<Atom, Charge>(frame);

        ASSERT_THAT(mass, NotNull()) << frame;
        ASSERT_THAT(charge, NotNull()) << frame;
        EXPECT_EQ(mass->size(), 3) << frame;
        EXPECT_EQ(charge->size(), 3) << frame;
    }
}

TEST_F(PropertyContainerTest, RemoveFrames)
{
    container.add<Atom, Mass>(true);

    for (size_t frame = 0; frame < 5; frame++)
    {
        container.add_frame();
    }

    Mass* mass_0 = container.get<Atom, Mass>(0);
    Mass* mass_2 = container.get<Atom, Mass>(2);
    Mass* mass_4 = container.get<Atom, Mass>(4);

    container.remove_frame(1);
    container.remove_frame(2); // Original frame 3

    Mass* new_mass_0 = container.get<Atom, Mass>(0);
    Mass* new_mass_1 = container.get<Atom, Mass>(1);
    Mass* new_mass_2 = container.get<Atom, Mass>(2);

    EXPECT_THAT(mass_0, NotNull());
    EXPECT_THAT(mass_2, NotNull());
    EXPECT_THAT(mass_4, NotNull());
    EXPECT_THAT(new_mass_0, NotNull());
    EXPECT_THAT(new_mass_1, NotNull());
    EXPECT_THAT(new_mass_2, NotNull());
    EXPECT_EQ(new_mass_0, mass_0);
    EXPECT_EQ(new_mass_1, mass_2);
    EXPECT_EQ(new_mass_2, mass_4);
}

TEST_F(PropertyContainerTest, RemoveInvalidFrame)
{
    container.add<Atom, Mass>(true);

    for (size_t frame = 0; frame < 3; frame++)
    {
        container.add_frame();
    }

    container.remove_frame(Frame());

    EXPECT_EQ(container.num_frames(), 3);
}

TEST_F(PropertyContainerTest, RemoveOutOfRangeFrame)
{
    container.add<Atom, Mass>(true);

    for (size_t frame = 0; frame < 3; frame++)
    {
        container.add_frame();
    }

    container.remove_frame(5);

    EXPECT_EQ(container.num_frames(), 3);
}
