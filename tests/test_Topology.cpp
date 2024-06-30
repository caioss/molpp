#include <molpp/internal/Topology.hpp>
#include <molpp/Common.hpp>

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

// Topology::link_entities should add a link between two entities
TEST_F(TopologyTest, link_entities)
{
    bool const result = topology.link_entities(atom1, residue1);

    EXPECT_TRUE(result);
    EXPECT_TRUE(topology.contains_link(atom1, residue1));
}

// Adding a link twice should return false and Topology should have only one link
TEST_F(TopologyTest, add_link_twice)
{
    ASSERT_TRUE(topology.link_entities(atom1, residue1));

    bool const result = topology.link_entities(atom1, residue1);

    EXPECT_FALSE(result);
    EXPECT_EQ(topology.count_links(atom1, MolecularEntityCategory::Residue), 1);
    EXPECT_EQ(topology.count_links(residue1, MolecularEntityCategory::Atom), 1);
}

// Topology::remove_link should remove a link between two entities
TEST_F(TopologyTest, remove_link)
{
    ASSERT_TRUE(topology.link_entities(atom1, residue1));

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

// Topology::find_link should return only the first link of a given entity
TEST_F(TopologyTest, find_link)
{
    ASSERT_TRUE(topology.link_entities(atom1, residue1));
    ASSERT_TRUE(topology.link_entities(atom2, residue1));

    std::optional<index_t> const result = topology.find_link(residue1, MolecularEntityCategory::Atom);

    ASSERT_TRUE(result.has_value());
    EXPECT_THAT(result.value(), AnyOf(Eq(atom1.index), Eq(atom2.index)));
}

// Topology::find_link should return std::nullopt if no link exists
TEST_F(TopologyTest, first_link_non_existent)
{
    ASSERT_TRUE(topology.link_entities(atom1, residue1));

    std::optional<index_t> const result = topology.find_link(atom2, MolecularEntityCategory::Residue);

    EXPECT_FALSE(result.has_value());
}

// Topology::all_links should return all links of a given entity
TEST_F(TopologyTest, all_links)
{
    ASSERT_TRUE(topology.link_entities(atom1, residue1));
    ASSERT_TRUE(topology.link_entities(atom2, residue1));

    std::vector<index_t> const result = topology.all_links(residue1, MolecularEntityCategory::Atom);

    EXPECT_THAT(result, UnorderedElementsAre(atom1.index, atom2.index));
}

// Topology::all_links should return an empty vector if no link exists
TEST_F(TopologyTest, all_links_non_existent)
{
    ASSERT_TRUE(topology.link_entities(atom1, residue1));

    std::vector<index_t> const result = topology.all_links(atom2, MolecularEntityCategory::Residue);

    EXPECT_TRUE(result.empty());
}

// Topology::count_links should return the number of links of a given entity
TEST_F(TopologyTest, count_links)
{
    ASSERT_TRUE(topology.link_entities(atom1, residue1));
    ASSERT_TRUE(topology.link_entities(atom2, residue1));

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
    ASSERT_TRUE(topology.link_entities(atom1, residue1));

    size_t const atom2_links = topology.count_links(atom2, MolecularEntityCategory::Residue);

    EXPECT_EQ(atom2_links, 0);
}

// Topology::contains_link should return true if a link exists
TEST_F(TopologyTest, contains_link)
{
    ASSERT_TRUE(topology.link_entities(atom1, residue1));

    bool const result = topology.contains_link(atom1, residue1);

    EXPECT_TRUE(result);
}

// Topology::contains_link should return false if no link exists
TEST_F(TopologyTest, contains_link_non_existent)
{
    ASSERT_TRUE(topology.link_entities(atom1, residue1));

    bool const result = topology.contains_link(atom2, residue1);

    EXPECT_FALSE(result);
}

// Topology::link_categories should return true when linking two categories
TEST_F(TopologyTest, link_categories)
{
    bool const result = topology.link_categories(MolecularEntityCategory::Atom, MolecularEntityCategory::Residue);

    EXPECT_TRUE(result);
}

// Linking two categories twice should return false
TEST_F(TopologyTest, link_categories_twice)
{
    ASSERT_TRUE(topology.link_categories(MolecularEntityCategory::Atom, MolecularEntityCategory::Residue));

    bool const result = topology.link_categories(MolecularEntityCategory::Atom, MolecularEntityCategory::Residue);

    EXPECT_FALSE(result);
}

// Topology::convert should convert a entity from one category to another one directly linked
TEST_F(TopologyTest, convert_directly_linked)
{
    ASSERT_TRUE(topology.link_categories(MolecularEntityCategory::Atom, MolecularEntityCategory::Residue));
    ASSERT_TRUE(topology.link_entities(atom1, residue1));
    ASSERT_TRUE(topology.link_entities(atom2, residue1));
    std::vector<index_t> const source_indices{residue1.index};

    std::vector<index_t> const result = topology.convert(source_indices, MolecularEntityCategory::Residue, MolecularEntityCategory::Atom);

    EXPECT_THAT(result, UnorderedElementsAre(atom1.index, atom2.index));
}

// Topology::convert should convert a entity from one category to another one not directly linked
TEST_F(TopologyTest, convert_indirectly_linked)
{
    Topology::MolecularEntityId const custom1{MolecularEntityCategory::Custom0, 2};
    Topology::MolecularEntityId const custom2{MolecularEntityCategory::Custom1, 1};

    ASSERT_TRUE(topology.link_categories(MolecularEntityCategory::Atom, MolecularEntityCategory::Residue));
    ASSERT_TRUE(topology.link_categories(MolecularEntityCategory::Residue, MolecularEntityCategory::Custom0));
    ASSERT_TRUE(topology.link_categories(MolecularEntityCategory::Custom0, MolecularEntityCategory::Custom1));
    ASSERT_TRUE(topology.link_entities(atom1, residue1));
    ASSERT_TRUE(topology.link_entities(atom2, residue1));
    ASSERT_TRUE(topology.link_entities(residue1, custom1));
    ASSERT_TRUE(topology.link_entities(custom1, custom2));
    std::vector<index_t> const source_indices{custom2.index};

    std::vector<index_t> const result = topology.convert(source_indices, custom2.category, MolecularEntityCategory::Atom);

    EXPECT_THAT(result, UnorderedElementsAre(atom1.index, atom2.index));
}

// Topology::convert should convert multiple entities from one category to another one not directly linked
TEST_F(TopologyTest, convert_multiple_indirectly_linked)
{
    Topology::MolecularEntityId const residue2{MolecularEntityCategory::Residue, 1};
    Topology::MolecularEntityId const custom1{MolecularEntityCategory::Custom0, 2};
    Topology::MolecularEntityId const custom2{MolecularEntityCategory::Custom0, 1};

    ASSERT_TRUE(topology.link_categories(MolecularEntityCategory::Atom, MolecularEntityCategory::Residue));
    ASSERT_TRUE(topology.link_categories(MolecularEntityCategory::Residue, MolecularEntityCategory::Custom0));
    ASSERT_TRUE(topology.link_entities(atom1, residue1));
    ASSERT_TRUE(topology.link_entities(atom2, residue2));
    ASSERT_TRUE(topology.link_entities(residue1, custom1));
    ASSERT_TRUE(topology.link_entities(residue2, custom2));
    std::vector<index_t> const source_indices{custom2.index, custom1.index};

    std::vector<index_t> const result = topology.convert(source_indices, custom1.category, MolecularEntityCategory::Atom);

    EXPECT_THAT(result, UnorderedElementsAre(atom1.index, atom2.index));
}

// Topology::convert should return an empty vector if there is a conversion path but no completly linked entities along the path
TEST_F(TopologyTest, convert_no_linked_entities)
{
    Topology::MolecularEntityId const custom1{MolecularEntityCategory::Custom0, 2};

    ASSERT_TRUE(topology.link_categories(MolecularEntityCategory::Atom, MolecularEntityCategory::Residue));
    ASSERT_TRUE(topology.link_categories(MolecularEntityCategory::Residue, MolecularEntityCategory::Custom0));
    ASSERT_TRUE(topology.link_entities(atom1, residue1));
    ASSERT_TRUE(topology.link_entities(atom2, residue1));
    std::vector<index_t> const source_indices{custom1.index};

    std::vector<index_t> const result = topology.convert(source_indices, custom1.category, MolecularEntityCategory::Atom);

    EXPECT_TRUE(result.empty());
}

// Topology::convert should return an empty vector if there is no conversion path
TEST_F(TopologyTest, convert_no_conversion_path)
{
    std::vector<index_t> const source_indices{atom1.index};

    std::vector<index_t> const result = topology.convert(source_indices, atom1.category, MolecularEntityCategory::Residue);

    EXPECT_TRUE(result.empty());
}

// Topology::convert should return an empty vector if the source category does not exist
TEST_F(TopologyTest, convert_non_existent_source_category)
{
    ASSERT_TRUE(topology.link_categories(MolecularEntityCategory::Atom, MolecularEntityCategory::Residue));
    ASSERT_TRUE(topology.link_entities(atom1, residue1));
    std::vector<index_t> const source_indices{0};

    std::vector<index_t> const result = topology.convert(source_indices, MolecularEntityCategory::Custom0, MolecularEntityCategory::Atom);

    EXPECT_TRUE(result.empty());
}

// Topology::convert should return an empty vector if the target category does not exist
TEST_F(TopologyTest, convert_non_existent_target_category)
{
    ASSERT_TRUE(topology.link_categories(MolecularEntityCategory::Atom, MolecularEntityCategory::Residue));
    ASSERT_TRUE(topology.link_entities(atom1, residue1));
    std::vector<index_t> const source_indices{residue1.index};

    std::vector<index_t> const result = topology.convert(source_indices, residue1.category, MolecularEntityCategory::Custom0);

    EXPECT_TRUE(result.empty());
}

// Topology::convert should return same indices if source and target category are the same
TEST_F(TopologyTest, convert_same_category)
{
    ASSERT_TRUE(topology.link_categories(MolecularEntityCategory::Atom, MolecularEntityCategory::Residue));
    ASSERT_TRUE(topology.link_entities(atom1, residue1));
    std::vector<index_t> const source_indices{residue1.index};

    std::vector<index_t> const result = topology.convert(source_indices, residue1.category, residue1.category);

    EXPECT_THAT(result, UnorderedElementsAre(residue1.index));
}
