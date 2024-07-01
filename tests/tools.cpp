#include "utils.hpp"
#include "tools/algorithms.hpp"
#include "tools/math.hpp"
#include "tools/SpatialSearch.hpp"
#include <molpp/Common.hpp>
#include <molpp/internal/SelIndices.hpp>
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <vector>
#include <numeric>

using namespace mol::internal;
using namespace testing;

TEST(DataStructures, SpatialSearch) {
    Eigen::Matrix3Xf points(3, 10);
    points << 0.394864, -1.92212, 0.507918, -1.10244, 2.86075, -2.68154, 2.28578, 2.88865, -1.65935, 2.21034, 0.665579, 1.90012, -0.467064, -2.63234, 2.86835, -1.38607, 0.375439, 2.77146, -1.7135, -1.00928, 0.0346084, -1.89917, -2.84799, -2.49805, 2.24333, -2.44597, 0.981083, -2.16332, -2.92808, 0.216672;

    SpatialSearch<Eigen::Matrix3Xf> search(points, 3.5);

    EXPECT_THAT(search.query(3, 3.0), UnorderedElementsAre(2, 3, 5, 8));
    EXPECT_THAT(search.query(3, 5.0), UnorderedElementsAre(0, 1, 2, 3, 5, 8, 9));

    EXPECT_THAT(search.pairs(3.0), UnorderedElementsAre(
        FieldsAre(3, 2, FloatNear(7.4041, 0.0001)),
        FieldsAre(5, 3, FloatNear(4.0494, 0.0001)),
        FieldsAre(6, 0, FloatNear(4.5555, 0.0001)),
        FieldsAre(6, 4, FloatNear(8.1384, 0.0001)),
        FieldsAre(8, 2, FloatNear(6.2570, 0.0001)),
        FieldsAre(8, 3, FloatNear(1.3393, 0.0001)),
        FieldsAre(8, 5, FloatNear(1.3845, 0.0001)),
        FieldsAre(9, 0, FloatNear(6.1342, 0.0001)),
        FieldsAre(9, 6, FloatNear(2.5074, 0.0001))));

    EXPECT_THAT(search.pairs(4.0), UnorderedElementsAre(
        FieldsAre(1, 0, FloatNear(10.6320, 0.0001)),
        FieldsAre(2, 0, FloatNear(9.6050, 0.0001)),
        FieldsAre(2, 1, FloatNear(12.4089, 0.0001)),
        FieldsAre(3, 2, FloatNear(7.4041, 0.0001)),
        FieldsAre(4, 0, FloatNear(15.8112, 0.0001)),
        FieldsAre(5, 1, FloatNear(11.6748, 0.0001)),
        FieldsAre(5, 2, FloatNear(11.1788, 0.0001)),
        FieldsAre(5, 3, FloatNear(4.0494, 0.0001)),
        FieldsAre(6, 0, FloatNear(4.5555, 0.0001)),
        FieldsAre(6, 4, FloatNear(8.1384, 0.0001)),
        FieldsAre(7, 0, FloatNear(15.4846, 0.0001)),
        FieldsAre(7, 6, FloatNear(15.9916, 0.0001)),
        FieldsAre(8, 1, FloatNear(14.1860, 0.0001)),
        FieldsAre(8, 2, FloatNear(6.2570, 0.0001)),
        FieldsAre(8, 3, FloatNear(1.3393, 0.0001)),
        FieldsAre(8, 5, FloatNear(1.3845, 0.0001)),
        FieldsAre(9, 0, FloatNear(6.1342, 0.0001)),
        FieldsAre(9, 2, FloatNear(12.5844, 0.0001)),
        FieldsAre(9, 6, FloatNear(2.5074, 0.0001))));
}

TEST(Math, Comparison) {
    EXPECT_TRUE(approximately_equal(95.1, 100.0, 0.05));
    EXPECT_FALSE(essentially_equal(95.1, 100.0, 0.05));

    EXPECT_TRUE(definitely_greater(106.0, 100.0, 0.05));
    EXPECT_FALSE(definitely_greater(95.1, 100.0, 0.05));

    EXPECT_FALSE(definitely_less(105.1, 100.0, 0.05));
    EXPECT_TRUE(definitely_less(94.9, 100.0, 0.05));
}
