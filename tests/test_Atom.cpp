#include "auxiliary.hpp"
#include "molpp/Atom.hpp"
#include "molpp/Residue.hpp"
#include "molpp/MolError.hpp"

#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace mol;
using namespace testing;

class AtomTest : public ::testing::Test
{
public:
    AtomTest()
    : data(create_moldata(3, 1, 1, 1, 1))
    , atom(1, 0, &data)
    , const_atom(1, 0, &data)
    , atom_no_frame(0, std::nullopt, &data)
    {}

    MolData data;
    Atom atom;
    Atom const const_atom;
    Atom atom_no_frame;
};

TEST_F(AtomTest, Comparison)
{
    EXPECT_TRUE(Atom(1, 0, &data) == Atom(1, 0, &data));
    EXPECT_FALSE(Atom(0, 0, &data) == Atom(1, 0, &data));
    EXPECT_FALSE(Atom(1, std::nullopt, &data) == Atom(1, 0, &data));
    EXPECT_FALSE(Atom(1, 0, &data) == Atom(1, 0, nullptr));
}

TEST_F(AtomTest, OperatorBool)
{
    EXPECT_FALSE(Atom());
    EXPECT_FALSE(Atom(4, 0, &data));
    EXPECT_FALSE(Atom(1, 0, nullptr));
    EXPECT_TRUE(atom_no_frame);
    EXPECT_TRUE(atom);
    EXPECT_TRUE(const_atom);
}

TEST_F(AtomTest, Frames)
{
    EXPECT_EQ(atom.frame(), 0);
    EXPECT_EQ(const_atom.frame(), 0);
    EXPECT_FALSE(atom_no_frame.frame());
}

TEST_F(AtomTest, SetValidFrame)
{
    atom.set_frame(0);
    EXPECT_EQ(atom.frame(), 0);
}

TEST_F(AtomTest, SetNullFrame)
{
    atom.set_frame(std::nullopt);
    EXPECT_FALSE(atom.frame());
}

TEST_F(AtomTest, SetInvalidFrame)
{
    EXPECT_THROW(atom.set_frame(1), MolError);
}

TEST_F(AtomTest, Residues)
{
    EXPECT_EQ(const_atom.resid(), 1);
    EXPECT_EQ(const_atom.residue_id(), 1);

    EXPECT_EQ(atom.resid(), 1);
    EXPECT_EQ(atom.residue_id(), 1);
    EXPECT_EQ(atom.residue(), Residue(1, 0, &data));
}

TEST_F(AtomTest, AtomicProperty)
{
    EXPECT_EQ(const_atom.atomic(), 1);
    EXPECT_EQ(atom.atomic(), 1);
    atom.set_atomic(2);
    EXPECT_EQ(atom.atomic(), 2);
}

TEST_F(AtomTest, OccupancyProperty)
{
    EXPECT_EQ(const_atom.occupancy(), 1.0);
    EXPECT_EQ(atom.occupancy(), 1.0);
    atom.set_occupancy(0.5);
    EXPECT_EQ(atom.occupancy(), 0.5);
}

TEST_F(AtomTest, TempfactorProperty)
{
    EXPECT_EQ(const_atom.tempfactor(), 1.0);
    EXPECT_EQ(atom.tempfactor(), 1.0);
    atom.set_tempfactor(0.5);
    EXPECT_EQ(atom.tempfactor(), 0.5);
}

TEST_F(AtomTest, MassProperty)
{
    EXPECT_EQ(const_atom.mass(), 1.0);
    EXPECT_EQ(atom.mass(), 1.0);
    atom.set_mass(0.5);
    EXPECT_EQ(atom.mass(), 0.5);
}

TEST_F(AtomTest, ChargeProperty)
{
    EXPECT_EQ(const_atom.charge(), 1.0);
    EXPECT_EQ(atom.charge(), 1.0);
    atom.set_charge(0.5);
    EXPECT_EQ(atom.charge(), 0.5);
}

TEST_F(AtomTest, RadiusProperty)
{
    EXPECT_EQ(const_atom.radius(), 1.0);
    EXPECT_EQ(atom.radius(), 1.0);
    atom.set_radius(0.5);
    EXPECT_EQ(atom.radius(), 0.5);
}

TEST_F(AtomTest, NameProperty)
{
    EXPECT_EQ(const_atom.name(), "B");
    EXPECT_EQ(atom.name(), "B");
    atom.set_name("CA");
    EXPECT_EQ(atom.name(), "CA");
}

TEST_F(AtomTest, TypeProperty)
{
    EXPECT_EQ(const_atom.type(), "B");
    EXPECT_EQ(atom.type(), "B");
    atom.set_type("C");
    EXPECT_EQ(atom.type(), "C");
}

TEST_F(AtomTest, AltLocProperty)
{
    EXPECT_EQ(const_atom.altloc(), "B");
    EXPECT_EQ(atom.altloc(), "B");
    atom.set_altloc("C");
    EXPECT_EQ(atom.altloc(), "C");
}

TEST_F(AtomTest, ResnameProperty)
{
    EXPECT_EQ(const_atom.resname(), "B");
    EXPECT_EQ(atom.resname(), "B");
}

TEST_F(AtomTest, SegidProperty)
{
    EXPECT_EQ(const_atom.segid(), "A");
    EXPECT_EQ(atom.segid(), "A");
}

TEST_F(AtomTest, ChainProperty)
{
    EXPECT_EQ(const_atom.chain(), "A");
    EXPECT_EQ(atom.chain(), "A");
}

TEST_F(AtomTest, Coordinates)
{
    EXPECT_THAT(const_atom.coords().reshaped(), ElementsAre(1, 1, 1));
    EXPECT_THAT(atom.coords().reshaped(), ElementsAre(1, 1, 1));
    atom.coords() *= 2;
    EXPECT_THAT(atom.coords().reshaped(), ElementsAre(2, 2, 2));
}

TEST_F(AtomTest, CoordinatesOnInvalidFrame)
{
    ASSERT_FALSE(atom_no_frame.frame());
    EXPECT_THROW(atom_no_frame.coords(), MolError);
}

TEST_F(AtomTest, AddValidBond)
{
    ASSERT_THAT(atom.add_bond(2), NotNull());
    ASSERT_THAT(atom.bond(2), NotNull());
    auto bond = atom.bond(2);
    ASSERT_THAT(bond, NotNull());
    EXPECT_EQ(bond->atom1(), 1);
    EXPECT_EQ(bond->atom2(), 2);
}

TEST_F(AtomTest, AddInvalidBond)
{
    EXPECT_THROW(atom.add_bond(1), MolError);
    EXPECT_THROW(atom.add_bond(3), MolError);
}

TEST_F(AtomTest, ReAddBond)
{
    atom.add_bond(2);
    EXPECT_EQ(atom.add_bond(2), atom.bond(2));
}

TEST_F(AtomTest, AddBondFromRValue)
{
    EXPECT_EQ(atom.add_bond(Atom(2, 0, &data)), atom.bond(Atom(2, 0, &data)));
}

TEST_F(AtomTest, BondsList)
{
    auto bonds_list = atom_no_frame.bonds();
    ASSERT_EQ(bonds_list.size(), 1);
    EXPECT_EQ(bonds_list[0]->atom1(), 0);
    EXPECT_EQ(bonds_list[0]->atom2(), 1);
}
