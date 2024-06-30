#include "utils.hpp"

#include "readers/StructureDetect.hpp"
#include <molpp/Common.hpp>
#include <molpp/internal/MolData.hpp>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace mol;
using namespace mol::internal;
using namespace testing;

//! Test fixture for StructureDetect class
class StructureDetectTest : public ::testing::Test
{
public:
    StructureDetectTest()
    : data{4}
    , topology{data.topology()}
    , detect{data}
    {}

    MolData data;
    Topology const& topology;
    StructureDetect detect;
};

// Register atom with same residue fields should render same residue
TEST_F(StructureDetectTest, register_same_residue)
{
    detect.register_atom(0, 1, "ALA", "AA", "A");
    detect.register_atom(1, 1, "ALA", "AA", "A");
    std::vector<index_t> const residue_0 = topology.all_links({EntityCategory::Residue, 0}, EntityCategory::Atom);

    EXPECT_THAT(residue_0, UnorderedElementsAre(0, 1));
}

// Register atom with different resid should render different residues
TEST_F(StructureDetectTest, register_different_resid)
{
    detect.register_atom(0, 1, "ALA", "AA", "A");
    detect.register_atom(1, 2, "ALA", "AA", "A");
    std::vector<index_t> const residue_0 = topology.all_links({EntityCategory::Residue, 0}, EntityCategory::Atom);
    std::vector<index_t> const residue_1 = topology.all_links({EntityCategory::Residue, 1}, EntityCategory::Atom);

    EXPECT_THAT(residue_0, UnorderedElementsAre(0));
    EXPECT_THAT(residue_1, UnorderedElementsAre(1));
}

// Register atom with different chain should render different residues and chains
TEST_F(StructureDetectTest, register_different_chain)
{
    detect.register_atom(0, 1, "ALA", "AA", "A");
    detect.register_atom(1, 1, "ALA", "AA", "B");
    std::vector<index_t> const residue_0 = topology.all_links({EntityCategory::Residue, 0}, EntityCategory::Atom);
    std::vector<index_t> const residue_1 = topology.all_links({EntityCategory::Residue, 1}, EntityCategory::Atom);
    std::vector<index_t> const chain_0 = topology.all_links({EntityCategory::Chain, 0}, EntityCategory::Residue);
    std::vector<index_t> const chain_1 = topology.all_links({EntityCategory::Chain, 1}, EntityCategory::Residue);

    EXPECT_THAT(residue_0, UnorderedElementsAre(0));
    EXPECT_THAT(residue_1, UnorderedElementsAre(1));
    EXPECT_THAT(chain_0, UnorderedElementsAre(0));
    EXPECT_THAT(chain_1, UnorderedElementsAre(1));
}

// Register atom with different resname should render different residues
TEST_F(StructureDetectTest, register_different_resname)
{
    detect.register_atom(0, 1, "ALA", "AA", "A");
    detect.register_atom(1, 1, "LYS", "AA", "A");
    std::vector<index_t> const residue_0 = topology.all_links({EntityCategory::Residue, 0}, EntityCategory::Atom);
    std::vector<index_t> const residue_1 = topology.all_links({EntityCategory::Residue, 1}, EntityCategory::Atom);

    EXPECT_THAT(residue_0, UnorderedElementsAre(0));
    EXPECT_THAT(residue_1, UnorderedElementsAre(1));
}

// Register atom with different segment should render different residues and segments
TEST_F(StructureDetectTest, register_different_segment)
{
    detect.register_atom(0, 1, "ALA", "AA", "A");
    detect.register_atom(1, 1, "ALA", "BB", "A");
    std::vector<index_t> const residue_0 = topology.all_links({EntityCategory::Residue, 0}, EntityCategory::Atom);
    std::vector<index_t> const residue_1 = topology.all_links({EntityCategory::Residue, 1}, EntityCategory::Atom);
    std::vector<index_t> const segment_0 = topology.all_links({EntityCategory::Segment, 0}, EntityCategory::Residue);
    std::vector<index_t> const segment_1 = topology.all_links({EntityCategory::Segment, 1}, EntityCategory::Residue);

    EXPECT_THAT(residue_0, UnorderedElementsAre(0));
    EXPECT_THAT(residue_1, UnorderedElementsAre(1));
    EXPECT_THAT(segment_0, UnorderedElementsAre(0));
    EXPECT_THAT(segment_1, UnorderedElementsAre(1));
}

// Register atom with non-empty chain should link chain and residues
TEST_F(StructureDetectTest, register_chain)
{
    detect.register_atom(0, 1, "ALA", "AA", "A");

    std::vector<index_t> const chain_0 = topology.all_links({EntityCategory::Chain, 0}, EntityCategory::Residue);

    EXPECT_THAT(chain_0, UnorderedElementsAre(0));
}

// Register same residue with different chains should render different residues linked to different chains
TEST_F(StructureDetectTest, register_same_residue_different_chain)
{
    detect.register_atom(0, 1, "ALA", "AA", "A");
    detect.register_atom(1, 1, "ALA", "AA", "B");

    std::vector<index_t> const chain_0 = topology.all_links({EntityCategory::Chain, 0}, EntityCategory::Residue);
    std::vector<index_t> const chain_1 = topology.all_links({EntityCategory::Chain, 1}, EntityCategory::Residue);

    EXPECT_THAT(chain_0, UnorderedElementsAre(0));
    EXPECT_THAT(chain_1, UnorderedElementsAre(1));
}

// Register atom with same chain should render same chain
TEST_F(StructureDetectTest, register_same_chain)
{
    detect.register_atom(0, 1, "ALA", "AA", "A");
    detect.register_atom(1, 1, "ALA", "AA", "A");

    std::vector<index_t> const chain_0 = topology.all_links({EntityCategory::Chain, 0}, EntityCategory::Residue);

    EXPECT_THAT(chain_0, UnorderedElementsAre(0));
}

// Register atom with empty chain should not link chain and residues
TEST_F(StructureDetectTest, register_empty_chain)
{
    detect.register_atom(0, 1, "ALA", "AA", "");

    std::vector<index_t> const chain_0 = topology.all_links({EntityCategory::Chain, 0}, EntityCategory::Residue);

    EXPECT_THAT(chain_0, IsEmpty());
}

// Register atom with non-empty segment should link segment and residues
TEST_F(StructureDetectTest, register_segment)
{
    detect.register_atom(0, 1, "ALA", "AA", "A");

    std::vector<index_t> const segment_0 = topology.all_links({EntityCategory::Segment, 0}, EntityCategory::Residue);

    EXPECT_THAT(segment_0, UnorderedElementsAre(0));
}

// Register same residue with different segments should render different residues linked to different segments
TEST_F(StructureDetectTest, register_same_residue_different_segment)
{
    detect.register_atom(0, 1, "ALA", "AA", "A");
    detect.register_atom(1, 1, "ALA", "BB", "A");

    std::vector<index_t> const segment_0 = topology.all_links({EntityCategory::Segment, 0}, EntityCategory::Residue);
    std::vector<index_t> const segment_1 = topology.all_links({EntityCategory::Segment, 1}, EntityCategory::Residue);

    EXPECT_THAT(segment_0, UnorderedElementsAre(0));
    EXPECT_THAT(segment_1, UnorderedElementsAre(1));
}

// Register atom with same segment should render same segment
TEST_F(StructureDetectTest, register_same_segment)
{
    detect.register_atom(0, 1, "ALA", "AA", "A");
    detect.register_atom(1, 1, "ALA", "AA", "A");

    std::vector<index_t> const segment_0 = topology.all_links({EntityCategory::Segment, 0}, EntityCategory::Residue);

    EXPECT_THAT(segment_0, UnorderedElementsAre(0));
}

// Register atom with empty segment should not link segment and residues
TEST_F(StructureDetectTest, register_empty_segment)
{
    detect.register_atom(0, 1, "ALA", "", "A");

    std::vector<index_t> const segment_0 = topology.all_links({EntityCategory::Segment, 0}, EntityCategory::Residue);

    EXPECT_THAT(segment_0, IsEmpty());
}

// Register atom with same residue fields out of order should re-use residues
TEST_F(StructureDetectTest, register_residues_out_of_order)
{
    detect.register_atom(0, 1, "ALA", "AA", "A");
    detect.register_atom(1, 1, "LYS", "AA", "A");
    detect.register_atom(2, 1, "PRO", "AA", "A");
    detect.register_atom(3, 1, "LYS", "AA", "A");
    detect.register_atom(4, 1, "PRO", "AA", "A");
    detect.register_atom(5, 1, "ALA", "AA", "A");

    std::vector<index_t> const residue_0 = topology.all_links({EntityCategory::Residue, 0}, EntityCategory::Atom);
    std::vector<index_t> const residue_1 = topology.all_links({EntityCategory::Residue, 1}, EntityCategory::Atom);
    std::vector<index_t> const residue_2 = topology.all_links({EntityCategory::Residue, 2}, EntityCategory::Atom);

    EXPECT_THAT(residue_0, UnorderedElementsAre(0, 5));
    EXPECT_THAT(residue_1, UnorderedElementsAre(1, 3));
    EXPECT_THAT(residue_2, UnorderedElementsAre(2, 4));
}

// Register atom with same chain out of order should re-use chains
TEST_F(StructureDetectTest, register_chains_out_of_order)
{
    detect.register_atom(0, 1, "ALA", "AA", "A");
    detect.register_atom(1, 1, "ALA", "AA", "B");
    detect.register_atom(2, 1, "ALA", "AA", "C");
    detect.register_atom(3, 1, "ALA", "AA", "B");
    detect.register_atom(4, 1, "ALA", "AA", "C");
    detect.register_atom(5, 1, "ALA", "AA", "A");

    std::vector<index_t> const chain_0 = topology.all_links({EntityCategory::Chain, 0}, EntityCategory::Residue);
    std::vector<index_t> const chain_1 = topology.all_links({EntityCategory::Chain, 1}, EntityCategory::Residue);
    std::vector<index_t> const chain_2 = topology.all_links({EntityCategory::Chain, 2}, EntityCategory::Residue);

    EXPECT_THAT(chain_0, UnorderedElementsAre(0));
    EXPECT_THAT(chain_1, UnorderedElementsAre(1));
    EXPECT_THAT(chain_2, UnorderedElementsAre(2));
}

// Register atom with same segment out of order should re-use segments
TEST_F(StructureDetectTest, register_segments_out_of_order)
{
    detect.register_atom(0, 1, "ALA", "AA", "A");
    detect.register_atom(1, 1, "ALA", "BB", "A");
    detect.register_atom(2, 1, "ALA", "CC", "A");
    detect.register_atom(3, 1, "ALA", "BB", "A");
    detect.register_atom(4, 1, "ALA", "CC", "A");
    detect.register_atom(5, 1, "ALA", "AA", "A");

    std::vector<index_t> const segment_0 = topology.all_links({EntityCategory::Segment, 0}, EntityCategory::Residue);
    std::vector<index_t> const segment_1 = topology.all_links({EntityCategory::Segment, 1}, EntityCategory::Residue);
    std::vector<index_t> const segment_2 = topology.all_links({EntityCategory::Segment, 2}, EntityCategory::Residue);

    EXPECT_THAT(segment_0, UnorderedElementsAre(0));
    EXPECT_THAT(segment_1, UnorderedElementsAre(1));
    EXPECT_THAT(segment_2, UnorderedElementsAre(2));
}

// StructureDetect::update_residue_data should set structures data
TEST_F(StructureDetectTest, update_residue_data)
{
    // Register atoms
    detect.register_atom(0, 1, "ALA", "BB", "A");
    detect.register_atom(1, 2, "LYS", "AA", "C");
    detect.register_atom(2, 2, "LYS", "AA", "C");
    detect.register_atom(3, 1, "ALA", "BB", "A");

    // Update data
    detect.update_residue_data(data);

    // Check residues data
    ResidueData const& residues_data = data.residues();
    EXPECT_EQ(residues_data.size(), 2);
    EXPECT_EQ(residues_data.id(0), 1);
    EXPECT_EQ(residues_data.id(1), 2);
    EXPECT_EQ(residues_data.name(0), "ALA");
    EXPECT_EQ(residues_data.name(1), "LYS");

    // Check chain data
    ChainData const& chain_data = data.chains();
    EXPECT_EQ(chain_data.size(), 2);
    EXPECT_EQ(chain_data.name(0), "A");
    EXPECT_EQ(chain_data.name(1), "C");

    // Check segment data
    SegmentData const& segment_data = data.segments();
    EXPECT_EQ(segment_data.size(), 2);
    EXPECT_EQ(segment_data.name(0), "BB");
    EXPECT_EQ(segment_data.name(1), "AA");
}

// StructureDetect::update_residue_data should not set empty chain data
TEST_F(StructureDetectTest, update_residue_data_empty_chain)
{
    // Register atoms
    detect.register_atom(0, 1, "ALA", "BB", "");
    detect.register_atom(1, 2, "LYS", "AA", "C");

    // Update data
    detect.update_residue_data(data);

    // Check chain data
    ChainData const& chain_data = data.chains();
    EXPECT_EQ(chain_data.size(), 1);
    EXPECT_EQ(chain_data.name(0), "C");
}

// StructureDetect::update_residue_data should not set empty segment data
TEST_F(StructureDetectTest, update_residue_data_empty_segment)
{
    // Register atoms
    detect.register_atom(0, 1, "ALA", "", "C");
    detect.register_atom(1, 2, "LYS", "AA", "C");

    // Update data
    detect.update_residue_data(data);

    // Check segment data
    SegmentData const& segment_data = data.segments();
    EXPECT_EQ(segment_data.size(), 1);
    EXPECT_EQ(segment_data.name(0), "AA");
}
