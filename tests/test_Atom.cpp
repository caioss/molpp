#include "auxiliary.hpp"
#include <molpp/Atom.hpp>
#include <molpp/Residue.hpp>
#include <molpp/MolError.hpp>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace mol;
using namespace testing;

//! Test fixture for Atom class.
class AtomTest : public ::testing::Test
{
public:
    AtomTest()
    : data(create_moldata(3, 1, 1, 1, 1))
    , atom(1, 0, data)
    , const_atom(1, 0, data)
    , null_frame_atom(0, std::nullopt, data)
    {}

    mol::internal::MolData data;
    Atom atom;
    Atom const const_atom;
    Atom null_frame_atom;
};

TEST_F(AtomTest, compare_atoms)
{
    EXPECT_TRUE(Atom(1, 0, data) == Atom(1, 0, data));
    EXPECT_FALSE(Atom(0, 0, data) == Atom(1, 0, data));
    EXPECT_FALSE(Atom(1, std::nullopt, data) == Atom(1, 0, data));
}

TEST_F(AtomTest, frames)
{
    EXPECT_EQ(atom.frame(), 0);
    EXPECT_EQ(const_atom.frame(), 0);
    EXPECT_FALSE(null_frame_atom.frame());
}

TEST_F(AtomTest, set_valid_frame)
{
    atom.set_frame(0);
    EXPECT_EQ(atom.frame(), 0);
}

TEST_F(AtomTest, set_null_frame)
{
    atom.set_frame(std::nullopt);
    EXPECT_FALSE(atom.frame());
}

TEST_F(AtomTest, set_invalid_frame)
{
    EXPECT_THROW(atom.set_frame(1), MolError);
}

TEST_F(AtomTest, residue_index)
{
    ASSERT_TRUE(atom.residue_index());
    EXPECT_EQ(atom.residue_index(), 1);
    EXPECT_EQ(const_atom.residue_index(), 1);
}

TEST_F(AtomTest, fetch_residue)
{
    ASSERT_TRUE(atom.residue());
    EXPECT_EQ(atom.residue(), Residue(1, 0, data));
}

// Atom::residue should return a null optional if the atom is not part of a residue.
TEST_F(AtomTest, fetch_null_residue)
{
    mol::internal::MolData data(1);
    mol::Atom atom(0, std::nullopt, data);

    EXPECT_FALSE(atom.residue());
}

TEST_F(AtomTest, change_residue)
{
    Residue new_res(0, 0, data);
    new_res.add_atom(atom);

    EXPECT_EQ(atom.residue_index(), 0);
}

TEST_F(AtomTest, atomic_property)
{
    EXPECT_EQ(const_atom.atomic_number(), 1);
    EXPECT_EQ(atom.atomic_number(), 1);
}

TEST_F(AtomTest, set_atomic_property)
{
    atom.set_atomic_number(2);

    EXPECT_EQ(atom.atomic_number(), 2);
}

TEST_F(AtomTest, occupancy_property)
{
    EXPECT_EQ(const_atom.occupancy(), 1.0);
    EXPECT_EQ(atom.occupancy(), 1.0);
}

TEST_F(AtomTest, set_occupancy_property)
{
    atom.set_occupancy(0.5);

    EXPECT_EQ(atom.occupancy(), 0.5);
}

TEST_F(AtomTest, temperature_factor_property)
{
    EXPECT_EQ(const_atom.temperature_factor(), 1.0);
    EXPECT_EQ(atom.temperature_factor(), 1.0);
}

TEST_F(AtomTest, set_temperature_factor_property)
{
    atom.set_temperature_factor(0.5);

    EXPECT_EQ(atom.temperature_factor(), 0.5);
}

TEST_F(AtomTest, mass_property)
{
    EXPECT_EQ(const_atom.mass(), 1.0);
    EXPECT_EQ(atom.mass(), 1.0);
}

TEST_F(AtomTest, set_mass_property)
{
    atom.set_mass(0.5);

    EXPECT_EQ(atom.mass(), 0.5);
}

TEST_F(AtomTest, charge_property)
{
    EXPECT_EQ(const_atom.charge(), 1.0);
    EXPECT_EQ(atom.charge(), 1.0);
}

TEST_F(AtomTest, set_charge_property)
{
    atom.set_charge(-0.5);

    EXPECT_EQ(atom.charge(), -0.5);
}

TEST_F(AtomTest, radius_property)
{
    EXPECT_EQ(const_atom.radius(), 1.0);
    EXPECT_EQ(atom.radius(), 1.0);
}

TEST_F(AtomTest, set_radius_property)
{
    atom.set_radius(0.5);

    EXPECT_EQ(atom.radius(), 0.5);
}

TEST_F(AtomTest, name_property)
{
    EXPECT_EQ(const_atom.name(), "B");
    EXPECT_EQ(atom.name(), "B");
}

TEST_F(AtomTest, set_name_property)
{
    atom.set_name("CA");

    EXPECT_EQ(atom.name(), "CA");
}

TEST_F(AtomTest, type_property)
{
    EXPECT_EQ(const_atom.type(), "B");
    EXPECT_EQ(atom.type(), "B");
}

TEST_F(AtomTest, set_type_property)
{
    atom.set_type("C");

    EXPECT_EQ(atom.type(), "C");
}

TEST_F(AtomTest, alternate_location_property)
{
    EXPECT_EQ(const_atom.alternate_location(), "B");
    EXPECT_EQ(atom.alternate_location(), "B");
}

TEST_F(AtomTest, set_alternate_location_property)
{
    atom.set_alternate_location("C");

    EXPECT_EQ(atom.alternate_location(), "C");
}

TEST_F(AtomTest, insertion_code_property)
{
    EXPECT_EQ(const_atom.insertion_code(), "B");
    EXPECT_EQ(atom.insertion_code(), "B");
}

TEST_F(AtomTest, set_insertion_code_property)
{
    atom.set_insertion_code("C");

    EXPECT_EQ(atom.insertion_code(), "C");
}

TEST_F(AtomTest, positions)
{
    EXPECT_THAT(const_atom.position().reshaped(), ElementsAre(1, 1, 1));
    EXPECT_THAT(atom.position().reshaped(), ElementsAre(1, 1, 1));
}

TEST_F(AtomTest, modify_positions)
{
    atom.position() *= 2;

    EXPECT_THAT(atom.position().reshaped(), ElementsAre(2, 2, 2));
}

TEST_F(AtomTest, positions_on_invalid_frame)
{
    ASSERT_FALSE(null_frame_atom.frame());
    EXPECT_THROW(null_frame_atom.position(), MolError);
}

TEST_F(AtomTest, add_valid_bond)
{
    auto added_bond = atom.add_bond(2);
    auto bond = atom.bond(2);

    ASSERT_THAT(added_bond, NotNull());
    ASSERT_THAT(bond, NotNull());
    EXPECT_EQ(added_bond, bond);
    EXPECT_EQ(added_bond->atom1(), 1);
    EXPECT_EQ(added_bond->atom2(), 2);
}

TEST_F(AtomTest, add_invalid_bond)
{
    EXPECT_THROW(atom.add_bond(1), MolError);
    EXPECT_THROW(atom.add_bond(3), MolError);
}

TEST_F(AtomTest, re_add_bond)
{
    atom.add_bond(2);

    EXPECT_EQ(atom.add_bond(2), atom.bond(2));
}

TEST_F(AtomTest, add_bond_from_r_value)
{
    EXPECT_EQ(atom.add_bond(Atom(2, 0, data)), atom.bond(Atom(2, 0, data)));
}

TEST_F(AtomTest, bonds_list)
{
    std::vector<std::pair<index_t, index_t>> bonds_indices;
    for (std::shared_ptr<Bond> bond : atom.bonds())
    {
        bonds_indices.push_back(std::make_pair<index_t, index_t>(bond->atom1(), bond->atom2()));
    }

    EXPECT_THAT(bonds_indices, UnorderedElementsAre(Pair(0, 1), Pair(1, 2)));
}

TEST_F(AtomTest, as_atom_indices)
{
    EXPECT_THAT(view2vector(atom.as_atom_indices()), ElementsAre(1));
    EXPECT_THAT(view2vector(const_atom.as_atom_indices()), ElementsAre(1));
    EXPECT_THAT(view2vector(null_frame_atom.as_atom_indices()), ElementsAre(0));
}
