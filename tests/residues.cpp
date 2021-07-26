#include "matchers.hpp"
#include "Atom.hpp"
#include "Residue.hpp"
#include "AtomSel.hpp"
#include "core/AtomData.hpp"
#include "core/ResidueData.hpp"
#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace mol;
using namespace mol::internal;
using namespace testing;

TEST(Residues, Residue) {
    // Construct data structures
    size_t const num_atoms { 3 };
    auto data = AtomData::create(num_atoms);
    ASSERT_THAT(data, NotNull());
    data->residues().resize(3);
    for (size_t i = 0; i < num_atoms; i++)
    {
        data->properties().residue(i) = i;
        data->residues().indices(i).insert(i);
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

    EXPECT_EQ(props.indices(0).size(), 0);
    EXPECT_EQ(props.resid(0), -1);
    EXPECT_EQ(props.resname(0), "");
    EXPECT_EQ(props.segid(0), "");
    EXPECT_EQ(props.chain(0), "");

    props.set(0, 1, "A", "B", "C");
    EXPECT_EQ(props.resid(0), 1);
    EXPECT_EQ(props.resname(0), "A");
    EXPECT_EQ(props.segid(0), "B");
    EXPECT_EQ(props.chain(0), "C");
}
