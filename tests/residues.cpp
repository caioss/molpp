#include "matchers.hpp"
#include "core/MolData.hpp"
#include "core/ResidueData.hpp"
#include "readers/ResidueDetect.hpp"
#include <molpp/Atom.hpp>
#include <molpp/Residue.hpp>
#include <molpp/AtomSel.hpp>
#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace mol;
using namespace mol::internal;
using namespace testing;

TEST(Residues, Residue) {
    // Construct data structures
    size_t const num_atoms { 3 };
    auto data = MolData::create(num_atoms);
    ASSERT_THAT(data, NotNull());
    data->residues().resize(3);
    for (size_t i = 0; i < num_atoms; i++)
    {
        data->properties().residue(i) = i;
        data->residues().add_atom(i, i);
    }

    // Comparison
    EXPECT_TRUE(Residue(1, 0, data) == Residue(1, 0, data));
    EXPECT_FALSE(Residue(0, 0, data) == Residue(1, 0, data));

    // Sizes
    ASSERT_EQ(data->size(), 3);

    /*
     * Properties
     */
    Residue res(1, 0, data);
    EXPECT_EQ(res.index(), 1);
    EXPECT_EQ(res.frame(), 0);

    res.set_resid(20);
    EXPECT_EQ(res.resid(), 20);

    res.set_resname("ARG");
    EXPECT_EQ(res.resname(), "ARG");

    res.set_segid("SEG1");
    EXPECT_EQ(res.segid(), "SEG1");

    res.set_chain("B");
    EXPECT_EQ(res.chain(), "B");

    /*
     * Atoms
     */
    ASSERT_EQ(res.size(), 1);
    auto sel = res.atoms();
    ASSERT_THAT(sel, NotNull());
    ASSERT_EQ(sel->size(), 1);
    EXPECT_EQ((*sel)[0].index(), 1);
    EXPECT_EQ((*sel)[0].residue_id(), 1);
    EXPECT_NE(Atom(0, 0, data).residue_id(), 1);
    EXPECT_NE(Atom(2, 0, data).residue_id(), 1);

    Atom atom = Atom(0, 0, data);
    res.add_atom(atom);
    res.add_atom(2);
    sel = res.atoms();
    ASSERT_THAT(sel, NotNull());
    ASSERT_EQ(sel->size(), 3);
    for (size_t i = 0; i < 3; i++)
    {
        EXPECT_EQ((*sel)[i].index(), i);
        EXPECT_EQ((*sel)[i].residue_id(), 1);
    }
}

TEST(Residues, ResidueData) {
    ResidueData props;
    EXPECT_EQ(props.size(), 0);
    props.resize(1);
    EXPECT_EQ(props.size(), 1);

    EXPECT_EQ(props.resid(0), -1);
    EXPECT_EQ(props.resname(0), "");
    EXPECT_EQ(props.segid(0), "");
    EXPECT_EQ(props.chain(0), "");

    props.set(0, 1, "A", "B", "C");
    EXPECT_EQ(props.resid(0), 1);
    EXPECT_EQ(props.resname(0), "A");
    EXPECT_EQ(props.segid(0), "B");
    EXPECT_EQ(props.chain(0), "C");

    EXPECT_EQ(props.size(0), 0);
    EXPECT_THAT(props.indices(0), UnorderedElementsAre());

    props.add_atom(0, 0);
    props.add_atom(0, 1);
    EXPECT_EQ(props.size(0), 2);
    EXPECT_THAT(props.indices(0), UnorderedElementsAre(0, 1));

    props.remove_atom(0, 1);
    EXPECT_EQ(props.size(0), 1);
    EXPECT_THAT(props.indices(0), UnorderedElementsAre(0));

    props.reset(0);
    EXPECT_EQ(props.size(0), 0);
    EXPECT_THAT(props.indices(0), UnorderedElementsAre());
}

TEST(Residues, ResidueDetect) {
    ResidueDetect detect;

    // New residues
    EXPECT_EQ(detect.register_atom(1, "ALA", "AA", "A"), 0);
    EXPECT_EQ(detect.register_atom(2, "ALA", "AA", "A"), 1);
    EXPECT_EQ(detect.register_atom(2, "MET", "AA", "A"), 2);
    EXPECT_EQ(detect.register_atom(2, "MET", "BB", "A"), 3);
    EXPECT_EQ(detect.register_atom(2, "MET", "BB", "B"), 4);

    // Repeated and non-sequential residues
    EXPECT_EQ(detect.register_atom(2, "ALA", "AA", "A"), 1);
    EXPECT_EQ(detect.register_atom(2, "MET", "BB", "B"), 4);
    EXPECT_EQ(detect.register_atom(2, "MET", "AA", "A"), 2);
    EXPECT_EQ(detect.register_atom(1, "ALA", "AA", "A"), 0);
    EXPECT_EQ(detect.register_atom(2, "MET", "BB", "A"), 3);

    // Expected data
    size_t index[10] = {0, 1, 2, 3, 4, 1, 4, 2, 0, 3};
    int resid[5] = {1, 2, 2, 2, 2};
    char resname[5][4] = {"ALA", "ALA", "MET", "MET", "MET"};
    char segid[5][3] = {"AA", "AA", "AA", "BB", "BB"};
    char chain[5][2] = {"A", "A", "A", "A", "B"};

    // Update MolData
    auto data = MolData::create(10);
    ASSERT_THAT(data, NotNull());
    for (size_t i = 0; i < 10; i++)
    {
        data->properties().residue(i) = index[i];
    }

    detect.update_residue_data(*data);
    ResidueData &residues_data = data->residues();

    // Check updates
    EXPECT_EQ(residues_data.size(), 5);

    for (size_t i = 0; i < 5; i++)
    {
        EXPECT_EQ(residues_data.resid(i), resid[i]) << "Residue " << i;
        EXPECT_EQ(residues_data.resname(i), resname[i]) << "Residue " << i;
        EXPECT_EQ(residues_data.segid(i), segid[i]) << "Residue " << i;
        EXPECT_EQ(residues_data.chain(i), chain[i]) << "Residue " << i;
    }

    EXPECT_THAT(residues_data.indices(0), UnorderedElementsAre(0, 8));
    EXPECT_THAT(residues_data.indices(1), UnorderedElementsAre(1, 5));
    EXPECT_THAT(residues_data.indices(2), UnorderedElementsAre(2, 7));
    EXPECT_THAT(residues_data.indices(3), UnorderedElementsAre(3, 9));
    EXPECT_THAT(residues_data.indices(4), UnorderedElementsAre(4, 6));
}
