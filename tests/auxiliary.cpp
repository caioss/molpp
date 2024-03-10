#include "auxiliary.hpp"
#include "matchers.hpp"
#include <molpp/Atom.hpp>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <optional>

using namespace testing;
using namespace mol;
using namespace mol::internal;

MolData create_moldata(size_t const num_res, size_t const num_res_atoms, size_t const num_chains, size_t const num_segments, size_t const num_frames)
{
    size_t const num_atoms { num_res * num_res_atoms };
    MolData data(num_atoms);
    data.add_entity<Atom>(num_atoms);
    AtomData& atom_data = data.atoms();
    ResidueData& res_data = data.residues();
    data.residues().resize(num_res);
    std::string const letters("ABCDEFGHIJKLMNOPQRSTUVWXYZ");

    // Properties
    Name* atom_name = data.add_property<Atom, Name>(false);
    Type* atom_type = data.add_property<Atom, Type>(false);
    AlternateLocation* atom_altloc = data.add_property<Atom, AlternateLocation>(false);
    InsertionCode* atom_insertion = data.add_property<Atom, InsertionCode>(false);
    AtomicNumber* atom_atomic = data.add_property<Atom, AtomicNumber>(false);
    Occupancy* atom_occupancy = data.add_property<Atom, Occupancy>(false);
    TemperatureFactor* atom_bfactor = data.add_property<Atom, TemperatureFactor>(false);
    Mass* atom_mass = data.add_property<Atom, Mass>(false);
    Charge* atom_charge = data.add_property<Atom, Charge>(false);
    Radius* atom_radius = data.add_property<Atom, Radius>(false);
    data.add_property<Atom, Position>(true);

    // Set atoms
    for (index_t atom_idx = 0; atom_idx < num_atoms; atom_idx++)
    {
        std::string const code = letters.substr(atom_idx % 26, 1);

        index_t const res_idx = atom_idx / num_res_atoms;
        atom_data.residue(atom_idx) = res_idx;
        res_data.add_atom(res_idx, atom_idx);

        // String properties
        atom_name->value(atom_idx) = code;
        atom_type->value(atom_idx) = code;
        atom_altloc->value(atom_idx) = code;
        atom_insertion->value(atom_idx) = code;

        // Numeric properties
        atom_atomic->value(atom_idx) = atom_idx;
        atom_occupancy->value(atom_idx) = atom_idx;
        atom_bfactor->value(atom_idx) = atom_idx;
        atom_mass->value(atom_idx) = atom_idx;
        atom_charge->value(atom_idx) = atom_idx;
        atom_radius->value(atom_idx) = atom_idx;
    }

    // Set residues
    for (index_t res_idx = 0; res_idx < num_res; res_idx++)
    {
        std::string const resname = letters.substr(res_idx % 26, 1);
        std::string const chain = letters.substr(res_idx % num_chains % 26, 1);
        std::string const segid = letters.substr(res_idx % num_segments % 26, 1);
        res_data.set(res_idx, res_idx, resname, segid, chain);
    }

    // Bonds between first atoms of consecutive residues
    for (index_t atom_idx = 0; atom_idx < num_atoms - num_res_atoms; atom_idx += num_res_atoms)
    {
        data.bonds().add_bond(atom_idx, atom_idx + num_res_atoms);
    }

    // Trajectory
    for (size_t frame_idx = 0; frame_idx < num_frames; frame_idx++)
    {
        Frame const frame = data.add_frame();
        Position* position_property = data.property_at<Atom, Position>(frame);
        Position::type& positions = position_property->positions();

        for (index_t atom_idx = 0; atom_idx < num_atoms; atom_idx++)
        {
            positions.col(atom_idx) << atom_idx, atom_idx, atom_idx;
        }
    }

    return data;
}

TEST(Auxiliary, create_moldata) {
    MolData data = create_moldata(3, 2, 2, 1, 2);
    ASSERT_EQ(data.entity_size<Atom>(), 6);

    // Atoms
    ASSERT_EQ(data.atoms().size(), 6);
    std::vector<mol::Atom> atoms;
    for (index_t i = 0; i < data.entity_size<Atom>(); ++i)
    {
        atoms.push_back(Atom(i, std::nullopt, &data));
    }

    EXPECT_THAT(atoms, Pointwise(Prop(&Atom::residue_id),
                                 {0ul, 0ul, 1ul, 1ul, 2ul, 2ul}));
    EXPECT_THAT(atoms, Pointwise(Prop(&Atom::resid),
                                 {0, 0, 1, 1, 2, 2}));

    // Residues
    EXPECT_EQ(data.residues().size(), 3);
}
