#include "utils.hpp"

#include <molpp/internal/MolData.hpp>
#include <molpp/ChainSel.hpp>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <array>
#include <utility>

using namespace testing;

//! Test fixture for ChainSel
class ChainSelTest : public ::testing::Test
{
public:
    ChainSelTest()
    : data(create_moldata(4, 1, 4, 4, 2))
    , selection{std::to_array<mol::index_t>({1, 2, 3}), data}
    , const_selection{std::to_array<mol::index_t>({1, 2, 3}), data}
    {
    }

    mol::internal::MolData data;
    mol::ChainSel selection;
    mol::ChainSel const const_selection;
};
