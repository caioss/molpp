#include <molpp/internal/Topology.hpp>
#include <molpp/MolppCore.hpp>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace mol;
using namespace mol::internal;
using namespace testing;

class TopologyTest : public testing::Test
{
public:
    TopologyTest()
    : atom1{MolecularEntityCategory::Atom, 0}
    , atom2{MolecularEntityCategory::Atom, 1}
    , residue1{MolecularEntityCategory::Residue, 0}
    {}

    Topology topology;
    Topology::MolecularEntityId atom1;
    Topology::MolecularEntityId atom2;
    Topology::MolecularEntityId residue1;
};

// Topology::add_link should add a link between two entities
TEST_F(TopologyTest, add_link)
{
    bool const result = topology.add_link(atom1, residue1);

    EXPECT_TRUE(result);
    EXPECT_TRUE(topology.contains_link(atom1, residue1));
}

// Adding a link twice should return false and Topology should have only one link
TEST_F(TopologyTest, add_link_twice)
{
    ASSERT_TRUE(topology.add_link(atom1, residue1));

    bool const result = topology.add_link(atom1, residue1);

    EXPECT_FALSE(result);
    EXPECT_EQ(topology.count_links(atom1, MolecularEntityCategory::Residue), 1);
    EXPECT_EQ(topology.count_links(residue1, MolecularEntityCategory::Atom), 1);
}

// Topology::remove_link should remove a link between two entities
TEST_F(TopologyTest, remove_link)
{
    ASSERT_TRUE(topology.add_link(atom1, residue1));

    bool const result = topology.remove_link(atom1, residue1);

    EXPECT_TRUE(result);
    EXPECT_FALSE(topology.contains_link(atom1, residue1));
}

// Removing a link that does not exist should return false
TEST_F(TopologyTest, remove_link_non_existent)
{
    bool const result = topology.remove_link(atom1, residue1);

    EXPECT_FALSE(result);
}

// Topology::first_link should return only the first link of a given entity
TEST_F(TopologyTest, first_link)
{
    ASSERT_TRUE(topology.add_link(atom1, residue1));
    ASSERT_TRUE(topology.add_link(atom2, residue1));

    std::optional<index_t> const result = topology.first_link(residue1, MolecularEntityCategory::Atom);

    ASSERT_TRUE(result.has_value());
    EXPECT_THAT(result.value(), AnyOf(Eq(atom1.index), Eq(atom2.index)));
}

// Topology::first_link should return std::nullopt if no link exists
TEST_F(TopologyTest, first_link_non_existent)
{
    ASSERT_TRUE(topology.add_link(atom1, residue1));

    std::optional<index_t> const result = topology.first_link(atom2, MolecularEntityCategory::Residue);

    EXPECT_FALSE(result.has_value());
}

// Topology::all_links should return all links of a given entity
TEST_F(TopologyTest, all_links)
{
    ASSERT_TRUE(topology.add_link(atom1, residue1));
    ASSERT_TRUE(topology.add_link(atom2, residue1));

    std::vector<index_t> const result = topology.all_links(residue1, MolecularEntityCategory::Atom);

    EXPECT_THAT(result, UnorderedElementsAre(atom1.index, atom2.index));
}

// Topology::all_links should return an empty vector if no link exists
TEST_F(TopologyTest, all_links_non_existent)
{
    ASSERT_TRUE(topology.add_link(atom1, residue1));

    std::vector<index_t> const result = topology.all_links(atom2, MolecularEntityCategory::Residue);

    EXPECT_TRUE(result.empty());
}

// Topology::count_links should return the number of links of a given entity
TEST_F(TopologyTest, count_links)
{
    ASSERT_TRUE(topology.add_link(atom1, residue1));
    ASSERT_TRUE(topology.add_link(atom2, residue1));

    size_t const atom1_links = topology.count_links(atom1, MolecularEntityCategory::Residue);
    size_t const atom2_links = topology.count_links(atom2, MolecularEntityCategory::Residue);
    size_t const residue1_links = topology.count_links(residue1, MolecularEntityCategory::Atom);

    EXPECT_EQ(atom1_links, 1);
    EXPECT_EQ(atom2_links, 1);
    EXPECT_EQ(residue1_links, 2);
}

// Topology::count_links should return 0 if no link exists
TEST_F(TopologyTest, count_links_non_existent)
{
    ASSERT_TRUE(topology.add_link(atom1, residue1));

    size_t const atom2_links = topology.count_links(atom2, MolecularEntityCategory::Residue);

    EXPECT_EQ(atom2_links, 0);
}

// Topology::contains_link should return true if a link exists
TEST_F(TopologyTest, contains_link)
{
    ASSERT_TRUE(topology.add_link(atom1, residue1));

    bool const result = topology.contains_link(atom1, residue1);

    EXPECT_TRUE(result);
}

// Topology::contains_link should return false if no link exists
TEST_F(TopologyTest, contains_link_non_existent)
{
    ASSERT_TRUE(topology.add_link(atom1, residue1));

    bool const result = topology.contains_link(atom2, residue1);

    EXPECT_FALSE(result);
}
