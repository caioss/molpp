#include "tools/SpatialSearch.hpp"

#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace testing;

class SpatialSearchTest : public testing::Test
{
public:
    SpatialSearchTest()
    : points(create_points())
    , search(points, 3.5)
    {}

    Eigen::Matrix3Xf points;
    mol::internal::SpatialSearch<Eigen::Matrix3Xf> search;

private:
    static Eigen::Matrix3Xf create_points()
    {
        Eigen::Matrix3Xf points(3, 10);
        points << 0.394864, -1.92212, 0.507918, -1.10244, 2.86075, -2.68154, 2.28578, 2.88865, -1.65935, 2.21034, 0.665579, 1.90012, -0.467064, -2.63234, 2.86835, -1.38607, 0.375439, 2.77146, -1.7135, -1.00928, 0.0346084, -1.89917, -2.84799, -2.49805, 2.24333, -2.44597, 0.981083, -2.16332, -2.92808, 0.216672;
        return points;
    }
};

// Calling SpatialSearch::query with a distance less than the cell size should return the indices of the points within the given distance
TEST_F(SpatialSearchTest, query_less_than_cell_size)
{
    EXPECT_THAT(search.query(3, 3.0), UnorderedElementsAre(2, 3, 5, 8));
}

// Calling SpatialSearch::query with a distance greater than the cell size should return the indices of the points within the given distance
TEST_F(SpatialSearchTest, query_greater_than_cell_size)
{
    EXPECT_THAT(search.query(3, 5.0), UnorderedElementsAre(0, 1, 2, 3, 5, 8, 9));
}

// Calling SpatialSearch::pairs with a distance less than the cell size should return all pairs of points within the given distance from each other and their squared distances
TEST_F(SpatialSearchTest, pairs_less_then_cell_size)
{
    EXPECT_THAT(search.pairs(3.0), UnorderedElementsAre(
        FieldsAre(3, 2, FloatEq(7.404130)),
        FieldsAre(5, 3, FloatEq(4.049458)),
        FieldsAre(6, 0, FloatEq(4.555558)),
        FieldsAre(6, 4, FloatEq(8.138463)),
        FieldsAre(8, 2, FloatEq(6.257067)),
        FieldsAre(8, 3, FloatEq(1.339341)),
        FieldsAre(8, 5, FloatEq(1.3845129)),
        FieldsAre(9, 0, FloatEq(6.134253)),
        FieldsAre(9, 6, FloatEq(2.507461))));
}

// Calling SpatialSearch::pairs with a distance greater than the cell size should return all pairs of points within the given distance from each other and their squared distances
TEST_F(SpatialSearchTest, pairs_greater_than_cell_size)
{
    EXPECT_THAT(search.pairs(4.0), UnorderedElementsAre(
        FieldsAre(1, 0, FloatEq(10.632004)),
        FieldsAre(2, 0, FloatEq(9.605034)),
        FieldsAre(2, 1, FloatEq(12.4089)),
        FieldsAre(3, 2, FloatEq(7.404130)),
        FieldsAre(4, 0, FloatEq(15.811245)),
        FieldsAre(5, 1, FloatEq(11.674753)),
        FieldsAre(5, 2, FloatEq(11.178834)),
        FieldsAre(5, 3, FloatEq(4.049458)),
        FieldsAre(6, 0, FloatEq(4.555558)),
        FieldsAre(6, 4, FloatEq(8.138463)),
        FieldsAre(7, 0, FloatEq(15.484592)),
        FieldsAre(7, 6, FloatEq(15.991639)),
        FieldsAre(8, 1, FloatEq(14.185953)),
        FieldsAre(8, 2, FloatEq(6.257067)),
        FieldsAre(8, 3, FloatEq(1.339341)),
        FieldsAre(8, 5, FloatEq(1.3845129)),
        FieldsAre(9, 0, FloatEq(6.134253)),
        FieldsAre(9, 2, FloatEq(12.584391)),
        FieldsAre(9, 6, FloatEq(2.507461))));
}
