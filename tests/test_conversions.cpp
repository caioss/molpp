#include "auxiliary.hpp"

#include <molpp/internal/MolData.hpp>
#include <molpp/Atom.hpp>
#include <molpp/Residue.hpp>
#include <molpp/Chain.hpp>
#include <molpp/Segment.hpp>
#include <molpp/AtomSel.hpp>
#include <molpp/ResidueSel.hpp>
#include <molpp/ChainSel.hpp>
#include <molpp/SegmentSel.hpp>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace testing;

//! Type traits class for type-parameterized tests
template<class SelType, class EntityType>
struct ConversionTraits
{
    using selection_type = SelType;
    using entity_type = EntityType;
};

//! Type-parameterized test fixture for convertions
template<class Traits>
class ConversionTest : public ::testing::Test
{
public:
    ConversionTest()
    : size(4)
    , data(create_moldata(size, 1, size, size, 1))
    {
    }

    size_t size;
    mol::internal::MolData data;
};

TYPED_TEST_SUITE_P(ConversionTest);

// Construct a selection from an entity
TYPED_TEST_P(ConversionTest, construct_from_entity)
{
    using Selection = typename TypeParam::selection_type;
    using Entity = typename TypeParam::entity_type;

    for (size_t index = 0; index < this->size; index++)
    {
        Selection selection{Entity{index, 0, this->data}};

        EXPECT_EQ(selection.size(), 1) << "Index " << index;
        EXPECT_THAT(selection.indices(), ElementsAre(index)) << "Index " << index;
    }
}

REGISTER_TYPED_TEST_SUITE_P(ConversionTest, construct_from_entity);

// All combinations
using ConversionTypes = ::testing::Types<
    ConversionTraits<mol::AtomSel, mol::Atom>,
    ConversionTraits<mol::AtomSel, mol::Residue>,
    ConversionTraits<mol::AtomSel, mol::Chain>,
    ConversionTraits<mol::AtomSel, mol::Segment>,
    ConversionTraits<mol::ResidueSel, mol::Atom>,
    ConversionTraits<mol::ResidueSel, mol::Residue>,
    ConversionTraits<mol::ResidueSel, mol::Chain>,
    ConversionTraits<mol::ResidueSel, mol::Segment>,
    ConversionTraits<mol::ChainSel, mol::Atom>,
    ConversionTraits<mol::ChainSel, mol::Residue>,
    ConversionTraits<mol::ChainSel, mol::Chain>,
    ConversionTraits<mol::ChainSel, mol::Segment>,
    ConversionTraits<mol::SegmentSel, mol::Atom>,
    ConversionTraits<mol::SegmentSel, mol::Residue>,
    ConversionTraits<mol::SegmentSel, mol::Chain>,
    ConversionTraits<mol::SegmentSel, mol::Segment>>;
INSTANTIATE_TYPED_TEST_SUITE_P(SelFromEntity, ConversionTest, ConversionTypes);
