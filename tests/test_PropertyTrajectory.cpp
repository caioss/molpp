#include <molpp/internal/PropertyTrajectory.hpp>
#include <molpp/Property.hpp>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace testing;
using namespace mol;
using namespace mol::internal;

struct MakePropertyMock
{
    MOCK_METHOD(void, call, (), ());
};

struct MakeProperty
{
    MakeProperty(NiceMock<MakePropertyMock>& mock)
    : mock{mock}
    {}

    std::unique_ptr<mol::Property> make()
    {
        mock.call();
        return std::make_unique<ContainerProperty<int>>();
    }

    NiceMock<MakePropertyMock>& mock;
};

class PropertyTrajectoryTest : public ::testing::Test
{
protected:
    PropertyTrajectoryTest()
    : make_property{make_property_mock}
    , static_trajectory{false, std::bind(&MakeProperty::make, make_property)}
    , dynamic_trajectory{true, std::bind(&MakeProperty::make, make_property)}
    {}

    NiceMock<MakePropertyMock> make_property_mock;
    MakeProperty make_property;
    PropertyTrajectory static_trajectory;
    PropertyTrajectory dynamic_trajectory;
};

TEST_F(PropertyTrajectoryTest, IsTimeBased)
{
    EXPECT_FALSE(static_trajectory.is_time_based());
    EXPECT_TRUE(dynamic_trajectory.is_time_based());
}

TEST_F(PropertyTrajectoryTest, DefaultFrames)
{
    EXPECT_EQ(static_trajectory.num_frames(), 0);
    EXPECT_THAT(static_trajectory.get(0), NotNull());

    EXPECT_EQ(dynamic_trajectory.num_frames(), 0);
    EXPECT_THAT(dynamic_trajectory.get(0), IsNull());
}

TEST_F(PropertyTrajectoryTest, DefaultSizeOnStatic)
{
    mol::Property* property = static_trajectory.get(0);
    ASSERT_THAT(property, NotNull());

    EXPECT_EQ(property->size(), 0);
}

TEST_F(PropertyTrajectoryTest, DefaultSizeOnDynamic)
{
    EXPECT_CALL(make_property_mock, call).Times(2);

    mol::Property* property0 = dynamic_trajectory.add_frame();
    mol::Property* property1 = dynamic_trajectory.add_frame();
    ASSERT_THAT(property0, NotNull());
    ASSERT_THAT(property1, NotNull());

    EXPECT_EQ(property0->size(), 0);
    EXPECT_EQ(property1->size(), 0);
}

TEST_F(PropertyTrajectoryTest, ResizeOnStatic)
{
    mol::Property* property = static_trajectory.get(0);
    ASSERT_THAT(property, NotNull());

    for (size_t i = 3; i < 6; i++)
    {
        static_trajectory.resize(i);

        EXPECT_EQ(property->size(), i);
    }
}

TEST_F(PropertyTrajectoryTest, ResizeOnDynamic)
{
    EXPECT_CALL(make_property_mock, call).Times(2);

    mol::Property* property0 = dynamic_trajectory.add_frame();
    mol::Property* property1 = dynamic_trajectory.add_frame();
    ASSERT_THAT(property0, NotNull());
    ASSERT_THAT(property1, NotNull());

    for (size_t i = 3; i < 6; i++)
    {
        dynamic_trajectory.resize(i);

        EXPECT_EQ(property0->size(), i);
        EXPECT_EQ(property1->size(), i);
    }
}

TEST_F(PropertyTrajectoryTest, ResizeSameSizeTwiceOnStatic)
{
    mol::Property* property = static_trajectory.get(0);
    ASSERT_THAT(property, NotNull());

    static_trajectory.resize(3);
    EXPECT_EQ(property->size(), 3);

    static_trajectory.resize(3);
    EXPECT_EQ(property->size(), 3);
}

TEST_F(PropertyTrajectoryTest, ResizeSameSizeTwiceOnDynamic)
{
    EXPECT_CALL(make_property_mock, call).Times(2);

    mol::Property* property0 = dynamic_trajectory.add_frame();
    mol::Property* property1 = dynamic_trajectory.add_frame();
    ASSERT_THAT(property0, NotNull());
    ASSERT_THAT(property1, NotNull());

    // First time
    dynamic_trajectory.resize(3);
    EXPECT_EQ(property0->size(), 3);
    EXPECT_EQ(property1->size(), 3);

    // Second time
    dynamic_trajectory.resize(3);
    EXPECT_EQ(property0->size(), 3);
    EXPECT_EQ(property1->size(), 3);
}

TEST_F(PropertyTrajectoryTest, AddFrameOnStatic)
{
    EXPECT_CALL(make_property_mock, call).Times(0);

    for (size_t i = 0; i < 3; i++)
    {
        mol::Property* property = static_trajectory.add_frame();

        EXPECT_THAT(property, IsNull());
        EXPECT_EQ(static_trajectory.num_frames(), 0);
    }
}

TEST_F(PropertyTrajectoryTest, AddFrameOnDynamic)
{
    EXPECT_CALL(make_property_mock, call).Times(3);

    for (size_t i = 0; i < 3; i++)
    {
        mol::Property* property = dynamic_trajectory.add_frame();

        EXPECT_THAT(property, NotNull());
        EXPECT_EQ(dynamic_trajectory.num_frames(), i + 1);
    }
}

TEST_F(PropertyTrajectoryTest, AddMultipleFramesAtOnceOnStatic)
{
    EXPECT_CALL(make_property_mock, call).Times(0);

    static_trajectory.add_frames(3);

    EXPECT_EQ(static_trajectory.num_frames(), 0);
}

TEST_F(PropertyTrajectoryTest, AddMultipleFramesOAtOncenDynamic)
{
    EXPECT_CALL(make_property_mock, call).Times(3);

    dynamic_trajectory.add_frames(3);

    EXPECT_EQ(dynamic_trajectory.num_frames(), 3);
}

TEST_F(PropertyTrajectoryTest, RemoveFrameOnStatic)
{
    EXPECT_CALL(make_property_mock, call).Times(0);

    mol::Property* property = static_trajectory.get(0);
    ASSERT_THAT(property, NotNull());

    static_trajectory.add_frames(3);
    static_trajectory.remove_frame(0);

    EXPECT_EQ(static_trajectory.get(0), property);
    EXPECT_EQ(static_trajectory.get(1), property);
    EXPECT_EQ(static_trajectory.get(2), property);
    EXPECT_EQ(static_trajectory.num_frames(), 0);
}

TEST_F(PropertyTrajectoryTest, RemoveFrameOnDynamic)
{
    EXPECT_CALL(make_property_mock, call).Times(3);

    mol::Property* property0 = dynamic_trajectory.add_frame();
    mol::Property* property1 = dynamic_trajectory.add_frame();
    mol::Property* property2 = dynamic_trajectory.add_frame();
    ASSERT_THAT(property0, NotNull());
    ASSERT_THAT(property1, NotNull());
    ASSERT_THAT(property2, NotNull());

    dynamic_trajectory.remove_frame(1);

    EXPECT_EQ(dynamic_trajectory.get(0), property0);
    EXPECT_EQ(dynamic_trajectory.get(1), property2);
    EXPECT_THAT(dynamic_trajectory.get(2), IsNull());
    EXPECT_EQ(dynamic_trajectory.num_frames(), 2);
}

TEST_F(PropertyTrajectoryTest, RemoveInvalidFrame)
{
    EXPECT_CALL(make_property_mock, call).Times(2);

    mol::Property* property0 = dynamic_trajectory.add_frame();
    mol::Property* property1 = dynamic_trajectory.add_frame();
    ASSERT_THAT(property0, NotNull());
    ASSERT_THAT(property1, NotNull());

    dynamic_trajectory.remove_frame(2);

    EXPECT_EQ(dynamic_trajectory.get(0), property0);
    EXPECT_EQ(dynamic_trajectory.get(1), property1);
    EXPECT_EQ(dynamic_trajectory.num_frames(), 2);
}

TEST_F(PropertyTrajectoryTest, GetOnStatic)
{
    EXPECT_CALL(make_property_mock, call).Times(0);

    mol::Property* property = static_trajectory.get(0);
    static_trajectory.add_frames(3);

    EXPECT_THAT(property, NotNull());
    EXPECT_EQ(static_trajectory.get(0), property);
    EXPECT_EQ(static_trajectory.get(1), property);
    EXPECT_EQ(static_trajectory.get(2), property);
    EXPECT_EQ(static_trajectory.get(3), property);
}

TEST_F(PropertyTrajectoryTest, GetOnDynamic)
{
    EXPECT_CALL(make_property_mock, call).Times(3);

    mol::Property* property0 = dynamic_trajectory.add_frame();
    mol::Property* property1 = dynamic_trajectory.add_frame();
    mol::Property* property2 = dynamic_trajectory.add_frame();
    ASSERT_THAT(property0, NotNull());
    ASSERT_THAT(property1, NotNull());
    ASSERT_THAT(property2, NotNull());

    EXPECT_EQ(dynamic_trajectory.get(0), property0);
    EXPECT_EQ(dynamic_trajectory.get(1), property1);
    EXPECT_EQ(dynamic_trajectory.get(2), property2);
    EXPECT_THAT(dynamic_trajectory.get(3), IsNull());
}
