#include "auxiliary.hpp"

#include <molpp/internal/MolData.hpp>
#include <molpp/SegmentSel.hpp>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <array>
#include <utility>

using namespace testing;

//! Test fixture for SegmentSel
class SegmentSelTest : public ::testing::Test
{
public:
    SegmentSelTest()
    : data(create_moldata(4, 1, 4, 4, 2))
    , selection{std::to_array<mol::index_t>({1, 2, 3}), data}
    , const_selection{std::to_array<mol::index_t>({1, 2, 3}), data}
    {
    }

    mol::internal::MolData data;
    mol::SegmentSel selection;
    mol::SegmentSel const const_selection;
};
