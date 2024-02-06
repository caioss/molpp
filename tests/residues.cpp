#include "matchers.hpp"
#include "auxiliary.hpp"
#include <molpp/internal/MolData.hpp>
#include <molpp/internal/ResidueData.hpp>
#include "readers/ResidueDetect.hpp"
#include <molpp/Atom.hpp>
#include <molpp/Residue.hpp>
#include <molpp/AtomSel.hpp>
#include <molpp/MolError.hpp>
#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace mol;
using namespace mol::internal;
using namespace testing;

TEST(Residues, Residue) {
    MolData data = create_moldata(3, 1, 1, 1, 1);

    // Comparison
    EXPECT_TRUE(Residue(1, 0, &data) == Residue(1, 0, &data));
    EXPECT_FALSE(Residue(0, 0, &data) == Residue(1, 0, &data));
    EXPECT_FALSE(Residue(1, std::nullopt, &data) == Residue(1, 0, &data));
    EXPECT_FALSE(Residue(1, 0, &data) == Residue(1, 0, nullptr));

    Residue res(1, 0, &data);
    Residue const const_res(1, 0, &data);

    /*
     * Constructors
     */
    EXPECT_FALSE(Residue());
    EXPECT_FALSE(Residue(4, 0, &data));
    EXPECT_FALSE(Residue(1, 0, nullptr));
    EXPECT_TRUE(res);
    EXPECT_TRUE(const_res);

    /*
     * Properties
     */
    EXPECT_EQ(res.index(), 1);
    EXPECT_EQ(res.frame(), 0);
    res.set_frame(std::nullopt);
    EXPECT_FALSE(res.frame());
    res.set_frame(0);
    EXPECT_EQ(res.frame(), 0);
    EXPECT_THROW(res.set_frame(1), MolError);

    res.set_resid(20);
    EXPECT_EQ(res.resid(), 20);

    res.set_resname("ARG");
    EXPECT_EQ(res.resname(), "ARG");

    res.set_segid("SEG1");
    EXPECT_EQ(res.segid(), "SEG1");

    res.set_chain("B");
    EXPECT_EQ(res.chain(), "B");

    /*
     * Bonds
     */
    auto bond_list = Residue(0, 0, &data).bonds();
    ASSERT_EQ(bond_list.size(), 1);
    EXPECT_EQ(bond_list[0]->atom1(), 0);
    EXPECT_EQ(bond_list[0]->atom2(), 1);

    /*
     * Addition/removal
     */
    Atom atom = Atom(0, 0, &data);
    res.add_atom(atom);
    res.add_atom(2);
    AtomSel atomsel(res);
    ASSERT_EQ(atomsel.size(), 3);
    for (index_t i = 0; i < 3; i++)
    {
        EXPECT_EQ(atomsel[i].index(), i);
        EXPECT_EQ(atomsel[i].residue_id(), 1);
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
    EXPECT_FALSE(props.indices(0));

    props.add_atom(0, 0);
    props.add_atom(0, 1);
    EXPECT_EQ(props.size(0), 2);
    EXPECT_THAT(view2vector(props.indices(0)), UnorderedElementsAre(0, 1));

    props.remove_atom(0, 1);
    EXPECT_EQ(props.size(0), 1);
    EXPECT_THAT(view2vector(props.indices(0)), UnorderedElementsAre(0));

    props.reset(0);
    EXPECT_EQ(props.size(0), 0);
    EXPECT_FALSE(props.indices(0));
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
    index_t index[10] = {0, 1, 2, 3, 4, 1, 4, 2, 0, 3};
    int resid[5] = {1, 2, 2, 2, 2};
    char resname[5][4] = {"ALA", "ALA", "MET", "MET", "MET"};
    char segid[5][3] = {"AA", "AA", "AA", "BB", "BB"};
    char chain[5][2] = {"A", "A", "A", "A", "B"};

    // Update MolData
    MolData data(10);
    data.add_entity<mol::Atom>(10);
    for (index_t i = 0; i < 10; i++)
    {
        data.atoms().residue(i) = index[i];
    }

    detect.update_residue_data(data);
    ResidueData &residues_data = data.residues();

    // Check updates
    EXPECT_EQ(residues_data.size(), 5);

    for (index_t i = 0; i < 5; i++)
    {
        EXPECT_EQ(residues_data.resid(i), resid[i]) << "Residue " << i;
        EXPECT_EQ(residues_data.resname(i), resname[i]) << "Residue " << i;
        EXPECT_EQ(residues_data.segid(i), segid[i]) << "Residue " << i;
        EXPECT_EQ(residues_data.chain(i), chain[i]) << "Residue " << i;
    }

    EXPECT_THAT(view2vector(residues_data.indices(0)), UnorderedElementsAre(0, 8));
    EXPECT_THAT(view2vector(residues_data.indices(1)), UnorderedElementsAre(1, 5));
    EXPECT_THAT(view2vector(residues_data.indices(2)), UnorderedElementsAre(2, 7));
    EXPECT_THAT(view2vector(residues_data.indices(3)), UnorderedElementsAre(3, 9));
    EXPECT_THAT(view2vector(residues_data.indices(4)), UnorderedElementsAre(4, 6));
}
