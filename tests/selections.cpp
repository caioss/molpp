#include "files.hpp"
#include "matchers.hpp"
#include <molpp/Atom.hpp>
#include <molpp/Residue.hpp>
#include <molpp/AtomSel.hpp>
#include <molpp/ResidueSel.hpp>
#include <molpp/MolError.hpp>
#include <molpp/internal/BaseSel.hpp>
#include <molpp/internal/Sel.hpp>
#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace mol;
using namespace mol::internal;
using namespace testing;

TEST(Selections, SelIndex) {
    /*
     * Construction
     */
    SelIndex all(5);
    EXPECT_EQ(all.size(), 5);
    EXPECT_THAT(all.indices(), ElementsAre(0, 1, 2, 3, 4));
    for (index_t i = 0; i < 5; ++i)
    {
        EXPECT_TRUE(all.contains(i)) << "index " << i;
    }
    EXPECT_FALSE(all.contains(6));
    EXPECT_TRUE(all.indices_begin() != all.indices_end());
    EXPECT_EQ(all.indices_end() - all.indices_begin(), 5);

    // Constructors accepting indexes
    std::vector<index_t> indices{4, 1, 1, 3};
    SelIndex some(indices, 5);
    SelIndex rvalue(std::vector<index_t>{4, 1, 1, 3}, 5);
    EXPECT_THROW(SelIndex(std::vector<index_t>{4, 1, 1, 3}, 1), MolError);
    EXPECT_EQ(some.size(), 3);
    EXPECT_EQ(rvalue.size(), 3);
    EXPECT_THAT(some.indices(), ElementsAre(1, 3, 4));
    EXPECT_THAT(rvalue.indices(), ElementsAre(1, 3, 4));
    for (index_t i : {1, 3, 4})
    {
        EXPECT_TRUE(some.contains(i)) << "Index " << i;
        EXPECT_TRUE(rvalue.contains(i)) << "Index " << i;
    }
    for (index_t i : {0, 2, 5})
    {
        EXPECT_FALSE(some.contains(i)) << "Index " << i;
        EXPECT_FALSE(rvalue.contains(i)) << "Index " << i;
    }
    EXPECT_TRUE(some.indices_begin() != some.indices_end());
    EXPECT_EQ(some.indices_end() - some.indices_begin(), 3);
    EXPECT_TRUE(rvalue.indices_begin() != rvalue.indices_end());
    EXPECT_EQ(rvalue.indices_end() - rvalue.indices_begin(), 3);
}

TEST(Selections, BaseSel) {
    // Data
    MolData* pdb_tiny = PDBFiles::tiny();
    MolData* pdb_traj = PDBFiles::traj();
    ASSERT_TRUE(pdb_tiny);
    ASSERT_TRUE(pdb_traj);

    /*
     * Construction
     */
    BaseSel all_sel(SelIndex(pdb_tiny->size()), pdb_tiny);
    EXPECT_EQ(all_sel.size(), pdb_tiny->size());
    EXPECT_FALSE(all_sel.frame());
    EXPECT_THAT(all_sel.indices(), ElementsAre(0, 1, 2, 3, 4, 5));
    for (index_t i = 0; i < 6; ++i)
    {
        EXPECT_TRUE(all_sel.contains(i)) << "index " << i;
    }
    EXPECT_FALSE(all_sel.contains(6));

    // Constructors accepting indexes
    std::vector<index_t> indices{4, 1, 1, 3};
    BaseSel some_sel(SelIndex(indices, pdb_tiny->size()), pdb_tiny);
    BaseSel rvalue_sel(SelIndex(std::vector<index_t>{4, 1, 1, 3}, pdb_tiny->size()), pdb_tiny);
    EXPECT_EQ(some_sel.size(), indices.size() - 1);
    EXPECT_EQ(rvalue_sel.size(), indices.size() - 1);
    EXPECT_FALSE(some_sel.frame());
    EXPECT_FALSE(rvalue_sel.frame());
    EXPECT_THAT(some_sel.indices(), ElementsAre(1, 3, 4));
    EXPECT_THAT(rvalue_sel.indices(), ElementsAre(1, 3, 4));
    for (index_t i : {1, 3, 4})
    {
        EXPECT_TRUE(some_sel.contains(i)) << "Index " << i;
        EXPECT_TRUE(rvalue_sel.contains(i)) << "Index " << i;
    }
    for (index_t i : {0, 2, 5})
    {
        EXPECT_FALSE(some_sel.contains(i)) << "Index " << i;
        EXPECT_FALSE(rvalue_sel.contains(i)) << "Index " << i;
    }

    /*
     * Trajectory
     */
    // Invalid Timestep
    EXPECT_THROW(all_sel.timestep(), MolError);

    BaseSel traj_sel(SelIndex(pdb_traj->size()), pdb_traj);
    ASSERT_TRUE(traj_sel.frame());
    EXPECT_EQ(traj_sel.frame(), 0);
    EXPECT_EQ(&(traj_sel.timestep()), &(pdb_traj->trajectory().timestep(0)));
    traj_sel.set_frame(3);
    ASSERT_TRUE(traj_sel.frame());
    EXPECT_EQ(traj_sel.frame(), 3);
    EXPECT_EQ(&(traj_sel.timestep()), &(pdb_traj->trajectory().timestep(3)));
    EXPECT_THROW(traj_sel.set_frame(4), MolError);
    ASSERT_TRUE(traj_sel.frame());
    EXPECT_EQ(traj_sel.frame(), 3);
    EXPECT_EQ(&(traj_sel.timestep()), &(pdb_traj->trajectory().timestep(3)));
}

TEST(Selections, Sel) {
    // Data
    MolData* pdb_tiny = PDBFiles::tiny();
    MolData* pdb_traj = PDBFiles::traj();
    ASSERT_TRUE(pdb_tiny);
    ASSERT_TRUE(pdb_traj);

    /*
     * Construction
     */
    Sel<Atom, AtomSel> all_sel(pdb_tiny);
    EXPECT_FALSE(all_sel.frame());
    EXPECT_EQ(all_sel.size(), pdb_tiny->size());
    EXPECT_THAT(all_sel.indices(), ElementsAre(0, 1, 2, 3, 4, 5));
    for (index_t i = 0; i < 6; ++i)
    {
        EXPECT_TRUE(all_sel.contains(i)) << "index " << i;
        EXPECT_EQ(all_sel.by_index(i).index(), i);
    }
    EXPECT_FALSE(all_sel.contains(6));
    EXPECT_THROW(all_sel.by_index(6), MolError);

    // Constructors accepting indexes
    std::vector<index_t> indices{4, 1, 1, 3};
    Sel<Atom, AtomSel> some_sel(indices, pdb_tiny);
    Sel<Atom, AtomSel> rvalue_sel(std::vector<index_t>{4, 1, 1, 3}, pdb_tiny);
    EXPECT_EQ(some_sel.size(), indices.size() - 1);
    EXPECT_EQ(rvalue_sel.size(), indices.size() - 1);
    EXPECT_FALSE(some_sel.frame());
    EXPECT_FALSE(rvalue_sel.frame());
    EXPECT_THAT(some_sel.indices(), ElementsAre(1, 3, 4));
    EXPECT_THAT(rvalue_sel.indices(), ElementsAre(1, 3, 4));
    for (index_t i : {1, 3, 4})
    {
        EXPECT_TRUE(some_sel.contains(i)) << "Index " << i;
        EXPECT_TRUE(rvalue_sel.contains(i)) << "Index " << i;
        EXPECT_EQ(some_sel.by_index(i).index(), i);
        EXPECT_EQ(rvalue_sel.by_index(i).index(), i);
    }
    for (index_t i : {0, 2, 5})
    {
        EXPECT_FALSE(some_sel.contains(i)) << "Index " << i;
        EXPECT_FALSE(rvalue_sel.contains(i)) << "Index " << i;
        EXPECT_THROW(some_sel.by_index(i), MolError);
        EXPECT_THROW(rvalue_sel.by_index(i), MolError);
    }

    /*
     * Conversions
     */
    Sel<Atom, AtomSel> from_atom(all_sel[0]);
    EXPECT_EQ(from_atom.frame(), all_sel.frame());
    EXPECT_THAT(from_atom.indices(), ElementsAre(0));

    Sel<Atom, AtomSel> from_residue(all_sel[0].residue());
    EXPECT_EQ(from_residue.frame(), all_sel.frame());
    EXPECT_THAT(from_residue.indices(), ElementsAre(0));

    Sel<Atom, AtomSel> from_atomsel{AtomSel(all_sel[0])};
    EXPECT_EQ(from_atomsel.frame(), all_sel.frame());
    EXPECT_THAT(from_atomsel.indices(), ElementsAre(0));

    Sel<Atom, AtomSel> from_residuesel{ResidueSel(all_sel[0])};
    EXPECT_EQ(from_residuesel.frame(), all_sel.frame());
    EXPECT_THAT(from_residuesel.indices(), ElementsAre(0));

    /*
     * Atoms getters
     */
    EXPECT_EQ(some_sel[0].index(), 1);
    EXPECT_FALSE(some_sel[0].frame());
    EXPECT_EQ(some_sel.at(0).index(), 1);
    EXPECT_FALSE(some_sel.at(0).frame());
    EXPECT_EQ(some_sel.by_index(1).index(), 1);
    EXPECT_EQ(rvalue_sel[0].index(), 1);
    EXPECT_EQ(some_sel[1].index(), 3);
    EXPECT_FALSE(some_sel[1].frame());
    EXPECT_EQ(some_sel.at(1).index(), 3);
    EXPECT_FALSE(some_sel.at(1).frame());
    EXPECT_EQ(some_sel.by_index(3).index(), 3);
    EXPECT_EQ(rvalue_sel[1].index(), 3);
    EXPECT_EQ(some_sel[2].index(), 4);
    EXPECT_FALSE(some_sel[2].frame());
    EXPECT_EQ(some_sel.at(2).index(), 4);
    EXPECT_FALSE(some_sel.at(2).frame());
    EXPECT_EQ(some_sel.by_index(4).index(), 4);
    EXPECT_EQ(rvalue_sel[2].index(), 4);
    EXPECT_THROW(some_sel.at(3).index(), MolError);
    EXPECT_THROW(some_sel.by_index(5), MolError);

    /*
     * Trajectory
     */
    // Invalid Timestep
    EXPECT_THROW(all_sel.timestep(), MolError);
    EXPECT_THROW(all_sel.coords(), MolError);

    Sel<Atom, AtomSel> traj_sel(pdb_traj);
    ASSERT_TRUE(traj_sel.frame());
    EXPECT_EQ(traj_sel.frame(), 0);
    EXPECT_EQ(&(traj_sel.timestep()), &(pdb_traj->trajectory().timestep(0)));
    EXPECT_THAT(traj_sel.coords().reshaped(), ElementsAre(1, -1, 0, 2, -2, 1));
    traj_sel.set_frame(3);
    ASSERT_TRUE(traj_sel.frame());
    EXPECT_EQ(traj_sel.frame(), 3);
    EXPECT_EQ(&(traj_sel.timestep()), &(pdb_traj->trajectory().timestep(3)));
    EXPECT_THROW(traj_sel.set_frame(4), MolError);
    ASSERT_TRUE(traj_sel.frame());
    EXPECT_EQ(traj_sel.frame(), 3);
    EXPECT_EQ(&(traj_sel.timestep()), &(pdb_traj->trajectory().timestep(3)));
    EXPECT_THAT(traj_sel.coords().reshaped(), ElementsAre(24, -24, 0, 48, -48, 24));
    traj_sel.coords().array() += 3;
    EXPECT_THAT(traj_sel.coords().reshaped(), ElementsAre(27, -21, 3, 51, -45, 27));
    // Undo changes
    traj_sel.coords().array() -= 3;

    // Correctly forwarding frames
    EXPECT_TRUE(traj_sel[0].frame());
    EXPECT_EQ(traj_sel[0].frame(), 3);

    /*
     * Iterators
     */
    std::vector<Atom> atoms(some_sel.begin(), some_sel.end());
    EXPECT_THAT(atoms, Pointwise(Prop(&Atom::name), {"NA", "O22", "ND"}));
    auto it = some_sel.begin();
    auto end = some_sel.end();
    EXPECT_EQ(end - it, 3);
    EXPECT_EQ((*it).name(), "NA");
    EXPECT_EQ((*++it).name(), "O22");
    it++;
    EXPECT_EQ((*it++).name(), "ND");
    EXPECT_TRUE(it == end);
    EXPECT_FALSE(it != end);

    // Correctly forwarding frames
    for (Atom atom : some_sel)
    {
        EXPECT_FALSE(atom.frame());
    }
    for (Atom atom : traj_sel)
    {
        EXPECT_TRUE(atom.frame());
    }

    /*
     * Bonds
     */
    EXPECT_THAT(all_sel.bonded().indices(), ElementsAre(0, 2, 3));
    EXPECT_EQ(all_sel.bonded().frame(), all_sel.frame());
    EXPECT_THAT(some_sel.bonded().indices(), ElementsAre(0, 3));
    EXPECT_EQ(some_sel.bonded().frame(), some_sel.frame());

    auto bonds = some_sel.bonds();
    EXPECT_EQ(bonds.size(), 1);
    EXPECT_EQ(bonds[0]->atom1(), 0);
    EXPECT_EQ(bonds[0]->atom2(), 3);

    bonds = all_sel.bonds();
    std::vector<index_t> bond_indices;
    bond_indices.reserve(4);
    for (auto b : bonds)
    {
        bond_indices.push_back(b->atom1());
        bond_indices.push_back(b->atom2());
    }
    EXPECT_THAT(bond_indices, UnorderedElementsAre(0, 2, 0, 3));
}

TEST(Selections, AtomSel) {
    // Data
    MolData* pdb_tiny = PDBFiles::tiny();
    MolData* pdb_traj = PDBFiles::traj();
    ASSERT_TRUE(pdb_tiny);
    ASSERT_TRUE(pdb_traj);

    /*
     * Construction
     */
    AtomSel all_sel(pdb_tiny);
    EXPECT_EQ(all_sel.size(), pdb_tiny->size());
    EXPECT_FALSE(all_sel.frame());
    EXPECT_THAT(all_sel.indices(), ElementsAre(0, 1, 2, 3, 4, 5));
    EXPECT_THAT(all_sel.atom_indices(), ElementsAre(0, 1, 2, 3, 4, 5));
    for (index_t i = 0; i < 6; ++i)
    {
        EXPECT_TRUE(all_sel.contains(i)) << "index " << i;
        EXPECT_EQ(all_sel.by_index(i).index(), i);
    }
    EXPECT_FALSE(all_sel.contains(6));
    EXPECT_THROW(all_sel.by_index(6), MolError);

    // Constructors accepting indexes
    std::vector<index_t> indices{4, 1, 1, 3};
    AtomSel some_sel(indices, pdb_tiny);
    AtomSel rvalue_sel(std::vector<index_t>{4, 1, 1, 3}, pdb_tiny);
    EXPECT_EQ(some_sel.size(), indices.size() - 1);
    EXPECT_EQ(rvalue_sel.size(), indices.size() - 1);
    EXPECT_FALSE(some_sel.frame());
    EXPECT_FALSE(rvalue_sel.frame());
    EXPECT_THAT(some_sel.indices(), ElementsAre(1, 3, 4));
    EXPECT_THAT(rvalue_sel.indices(), ElementsAre(1, 3, 4));
    EXPECT_THAT(some_sel.atom_indices(), ElementsAre(1, 3, 4));
    EXPECT_THAT(rvalue_sel.atom_indices(), ElementsAre(1, 3, 4));
    for (index_t i : {1, 3, 4})
    {
        EXPECT_TRUE(some_sel.contains(i)) << "Index " << i;
        EXPECT_TRUE(rvalue_sel.contains(i)) << "Index " << i;
        EXPECT_EQ(some_sel.by_index(i).index(), i);
        EXPECT_EQ(rvalue_sel.by_index(i).index(), i);
    }
    for (index_t i : {0, 2, 5})
    {
        EXPECT_FALSE(some_sel.contains(i)) << "Index " << i;
        EXPECT_FALSE(rvalue_sel.contains(i)) << "Index " << i;
        EXPECT_THROW(some_sel.by_index(i), MolError);
        EXPECT_THROW(rvalue_sel.by_index(i), MolError);
    }

    /*
     * Conversions
     */
    AtomSel from_atom(all_sel[0]);
    EXPECT_EQ(from_atom.frame(), all_sel.frame());
    EXPECT_THAT(from_atom.indices(), ElementsAre(0));

    AtomSel from_residue(all_sel[0].residue());
    EXPECT_EQ(from_residue.frame(), all_sel.frame());
    EXPECT_THAT(from_residue.indices(), ElementsAre(0));

    AtomSel from_atomsel(all_sel);
    EXPECT_EQ(from_atomsel.frame(), all_sel.frame());
    EXPECT_THAT(from_atomsel.indices(), ContainerEq(all_sel.indices()));

    AtomSel from_residuesel{ResidueSel(all_sel[0])};
    EXPECT_EQ(from_residuesel.frame(), all_sel.frame());
    EXPECT_THAT(from_residuesel.indices(), ElementsAre(0));

    /*
     * Atoms getters
     */
    EXPECT_EQ(some_sel[0].index(), 1);
    EXPECT_FALSE(some_sel[0].frame());
    EXPECT_EQ(some_sel.at(0).index(), 1);
    EXPECT_FALSE(some_sel.at(0).frame());
    EXPECT_EQ(some_sel.by_index(1).index(), 1);
    EXPECT_EQ(rvalue_sel[0].index(), 1);
    EXPECT_EQ(some_sel[1].index(), 3);
    EXPECT_FALSE(some_sel[1].frame());
    EXPECT_EQ(some_sel.at(1).index(), 3);
    EXPECT_FALSE(some_sel.at(1).frame());
    EXPECT_EQ(some_sel.by_index(3).index(), 3);
    EXPECT_EQ(rvalue_sel[1].index(), 3);
    EXPECT_EQ(some_sel[2].index(), 4);
    EXPECT_FALSE(some_sel[2].frame());
    EXPECT_EQ(some_sel.at(2).index(), 4);
    EXPECT_FALSE(some_sel.at(2).frame());
    EXPECT_EQ(some_sel.by_index(4).index(), 4);
    EXPECT_EQ(rvalue_sel[2].index(), 4);
    EXPECT_THROW(some_sel.at(3).index(), MolError);
    EXPECT_THROW(some_sel.by_index(5), MolError);

    /*
     * Trajectory
     */
    // Invalid Timestep
    EXPECT_THROW(all_sel.timestep(), MolError);
    EXPECT_THROW(all_sel.coords(), MolError);

    AtomSel traj_sel(pdb_traj);
    ASSERT_TRUE(traj_sel.frame());
    EXPECT_EQ(traj_sel.frame(), 0);
    EXPECT_EQ(&(traj_sel.timestep()), &(pdb_traj->trajectory().timestep(0)));
    EXPECT_THAT(traj_sel.coords().reshaped(), ElementsAre(1, -1, 0, 2, -2, 1));
    traj_sel.set_frame(3);
    ASSERT_TRUE(traj_sel.frame());
    EXPECT_EQ(traj_sel.frame(), 3);
    EXPECT_EQ(&(traj_sel.timestep()), &(pdb_traj->trajectory().timestep(3)));
    EXPECT_THROW(traj_sel.set_frame(4), MolError);
    ASSERT_TRUE(traj_sel.frame());
    EXPECT_EQ(traj_sel.frame(), 3);
    EXPECT_EQ(&(traj_sel.timestep()), &(pdb_traj->trajectory().timestep(3)));
    EXPECT_THAT(traj_sel.coords().reshaped(), ElementsAre(24, -24, 0, 48, -48, 24));
    traj_sel.coords().array() += 3;
    EXPECT_THAT(traj_sel.coords().reshaped(), ElementsAre(27, -21, 3, 51, -45, 27));
    // Undo changes
    traj_sel.coords().array() -= 3;

    // Correctly forwarding frames
    EXPECT_TRUE(traj_sel[0].frame());
    EXPECT_EQ(traj_sel[0].frame(), 3);

    /*
     * Iterators
     */
    std::vector<Atom> atoms(some_sel.begin(), some_sel.end());
    EXPECT_THAT(atoms, Pointwise(Prop(&Atom::name), {"NA", "O22", "ND"}));
    auto it = some_sel.begin();
    auto end = some_sel.end();
    EXPECT_EQ(end - it, 3);
    EXPECT_EQ((*it).name(), "NA");
    EXPECT_EQ((*++it).name(), "O22");
    it++;
    EXPECT_EQ((*it++).name(), "ND");
    EXPECT_TRUE(it == end);
    EXPECT_FALSE(it != end);

    // Correctly forwarding frames
    for (Atom atom : some_sel)
    {
        EXPECT_FALSE(atom.frame());
    }
    for (Atom atom : traj_sel)
    {
        EXPECT_TRUE(atom.frame());
    }

    /*
     * Bonds
     */
    EXPECT_THAT(all_sel.bonded().indices(), ElementsAre(0, 2, 3));
    EXPECT_EQ(all_sel.bonded().frame(), all_sel.frame());
    EXPECT_THAT(some_sel.bonded().indices(), ElementsAre(0, 3));
    EXPECT_EQ(some_sel.bonded().frame(), some_sel.frame());

    auto bonds = some_sel.bonds();
    EXPECT_EQ(bonds.size(), 1);
    EXPECT_EQ(bonds[0]->atom1(), 0);
    EXPECT_EQ(bonds[0]->atom2(), 3);

    bonds = all_sel.bonds();
    std::vector<index_t> bond_indices;
    bond_indices.reserve(4);
    for (auto b : bonds)
    {
        bond_indices.push_back(b->atom1());
        bond_indices.push_back(b->atom2());
    }
    EXPECT_THAT(bond_indices, UnorderedElementsAre(0, 2, 0, 3));
}

TEST(Selections, ResidueSel) {
    // Data
    MolData* pdb_tiny = PDBFiles::tiny();
    MolData* pdb_traj = PDBFiles::traj();
    ASSERT_TRUE(pdb_tiny);
    ASSERT_TRUE(pdb_traj);

    /*
     * Construction
     */
    ResidueSel all_sel(pdb_tiny);
    EXPECT_EQ(all_sel.size(), 5);
    EXPECT_FALSE(all_sel.frame());
    EXPECT_THAT(all_sel.indices(), ElementsAre(0, 1, 2, 3, 4));
    EXPECT_THAT(all_sel.atom_indices(), UnorderedElementsAre(0, 1, 2, 3, 4, 5));
    for (index_t i = 0; i < 5; ++i)
    {
        EXPECT_TRUE(all_sel.contains(i)) << "index " << i;
        EXPECT_EQ(all_sel.by_index(i).index(), i);
    }
    EXPECT_FALSE(all_sel.contains(6));
    EXPECT_THROW(all_sel.by_index(6), MolError);

    // Constructors accepting indexes
    std::vector<index_t> indices{4, 1, 1, 3};
    ResidueSel some_sel(indices, pdb_tiny);
    ResidueSel rvalue_sel(std::vector<index_t>{4, 1, 1, 3}, pdb_tiny);
    EXPECT_EQ(some_sel.size(), indices.size() - 1);
    EXPECT_EQ(rvalue_sel.size(), indices.size() - 1);
    EXPECT_FALSE(some_sel.frame());
    EXPECT_FALSE(rvalue_sel.frame());
    EXPECT_THAT(some_sel.indices(), ElementsAre(1, 3, 4));
    EXPECT_THAT(rvalue_sel.indices(), ElementsAre(1, 3, 4));
    EXPECT_THAT(some_sel.atom_indices(), UnorderedElementsAre(1, 3, 4, 5));
    EXPECT_THAT(rvalue_sel.atom_indices(), UnorderedElementsAre(1, 3, 4, 5));
    for (index_t i : {1, 3, 4})
    {
        EXPECT_TRUE(some_sel.contains(i)) << "Index " << i;
        EXPECT_TRUE(rvalue_sel.contains(i)) << "Index " << i;
        EXPECT_EQ(some_sel.by_index(i).index(), i);
        EXPECT_EQ(rvalue_sel.by_index(i).index(), i);
    }
    for (index_t i : {0, 2, 5})
    {
        EXPECT_FALSE(some_sel.contains(i)) << "Index " << i;
        EXPECT_FALSE(rvalue_sel.contains(i)) << "Index " << i;
        EXPECT_THROW(some_sel.by_index(i), MolError);
        EXPECT_THROW(rvalue_sel.by_index(i), MolError);
    }

    /*
     * Conversions
     */
    ResidueSel from_atomsel{AtomSel(all_sel)};
    EXPECT_EQ(from_atomsel.frame(), all_sel.frame());
    EXPECT_THAT(from_atomsel.indices(), ElementsAre(0, 1, 2, 3, 4));

    ResidueSel from_residuesel(all_sel);
    EXPECT_EQ(from_residuesel.frame(), all_sel.frame());
    EXPECT_THAT(from_residuesel.indices(), ContainerEq(all_sel.indices()));

    ResidueSel from_atom(from_atomsel[0]);
    EXPECT_EQ(from_atom.frame(), all_sel.frame());
    EXPECT_THAT(from_atom.indices(), ElementsAre(0));

    ResidueSel from_residue(all_sel[0]);
    EXPECT_EQ(from_residue.frame(), all_sel.frame());
    EXPECT_THAT(from_residue.indices(), ElementsAre(0));

    /*
     * Indexing
     */
    EXPECT_EQ(some_sel[0].index(), 1);
    EXPECT_FALSE(some_sel[0].frame());
    EXPECT_EQ(some_sel.at(0).index(), 1);
    EXPECT_FALSE(some_sel.at(0).frame());
    EXPECT_EQ(some_sel.by_index(1).index(), 1);
    EXPECT_EQ(rvalue_sel[0].index(), 1);
    EXPECT_EQ(some_sel[1].index(), 3);
    EXPECT_FALSE(some_sel[1].frame());
    EXPECT_EQ(some_sel.at(1).index(), 3);
    EXPECT_FALSE(some_sel.at(1).frame());
    EXPECT_EQ(some_sel.by_index(3).index(), 3);
    EXPECT_EQ(rvalue_sel[1].index(), 3);
    EXPECT_EQ(some_sel[2].index(), 4);
    EXPECT_FALSE(some_sel[2].frame());
    EXPECT_EQ(some_sel.at(2).index(), 4);
    EXPECT_FALSE(some_sel.at(2).frame());
    EXPECT_EQ(some_sel.by_index(4).index(), 4);
    EXPECT_EQ(rvalue_sel[2].index(), 4);
    EXPECT_THROW(some_sel.at(3).index(), MolError);
    EXPECT_THROW(some_sel.by_index(5), MolError);

    /*
     * Trajectory
     */
    // Invalid Timestep
    EXPECT_THROW(all_sel.timestep(), MolError);
    EXPECT_THROW(all_sel.coords(), MolError);

    ResidueSel traj_sel(pdb_traj);
    ASSERT_TRUE(traj_sel.frame());
    EXPECT_EQ(traj_sel.frame(), 0);
    EXPECT_EQ(&(traj_sel.timestep()), &(pdb_traj->trajectory().timestep(0)));
    EXPECT_THAT(traj_sel.coords().reshaped(), ElementsAre(1, -1, 0, 2, -2, 1));
    traj_sel.set_frame(3);
    ASSERT_TRUE(traj_sel.frame());
    EXPECT_EQ(traj_sel.frame(), 3);
    EXPECT_EQ(&(traj_sel.timestep()), &(pdb_traj->trajectory().timestep(3)));
    EXPECT_THROW(traj_sel.set_frame(4), MolError);
    ASSERT_TRUE(traj_sel.frame());
    EXPECT_EQ(traj_sel.frame(), 3);
    EXPECT_EQ(&(traj_sel.timestep()), &(pdb_traj->trajectory().timestep(3)));
    EXPECT_THAT(traj_sel.coords().reshaped(), ElementsAre(24, -24, 0, 48, -48, 24));
    traj_sel.coords().array() += 3;
    EXPECT_THAT(traj_sel.coords().reshaped(), ElementsAre(27, -21, 3, 51, -45, 27));
    // Undo changes
    traj_sel.coords().array() -= 3;

    // Correctly forwarding frames
    EXPECT_TRUE(traj_sel[0].frame());
    EXPECT_EQ(traj_sel[0].frame(), 3);

    /*
     * Iterators
     */
    std::vector<Residue> residues(some_sel.begin(), some_sel.end());
    EXPECT_THAT(residues, Pointwise(Prop(&Residue::resid), {339, 801, 85}));
    auto it = some_sel.begin();
    auto end = some_sel.end();
    EXPECT_EQ(end - it, 3);
    EXPECT_EQ((*it).resid(), 339);
    EXPECT_EQ((*++it).resid(), 801);
    it++;
    EXPECT_EQ((*it++).resid(), 85);
    EXPECT_TRUE(it == end);
    EXPECT_FALSE(it != end);

    // Correctly forwarding frames
    for (Residue res : some_sel)
    {
        EXPECT_FALSE(res.frame());
    }
    for (Residue res : traj_sel)
    {
        EXPECT_TRUE(res.frame());
    }

    /*
     * Bonds
     */
    EXPECT_THAT(all_sel.bonded().indices(), ElementsAre(0, 2, 3));
    EXPECT_EQ(all_sel.bonded().frame(), all_sel.frame());
    EXPECT_THAT(some_sel.bonded().indices(), ElementsAre(0, 3));
    EXPECT_EQ(some_sel.bonded().frame(), some_sel.frame());

    auto bonds = some_sel.bonds();
    EXPECT_EQ(bonds.size(), 1);
    EXPECT_EQ(bonds[0]->atom1(), 0);
    EXPECT_EQ(bonds[0]->atom2(), 3);

    bonds = all_sel.bonds();
    std::vector<index_t> bond_indices;
    bond_indices.reserve(4);
    for (auto b : bonds)
    {
        bond_indices.push_back(b->atom1());
        bond_indices.push_back(b->atom2());
    }
    EXPECT_THAT(bond_indices, UnorderedElementsAre(0, 2, 0, 3));
}
