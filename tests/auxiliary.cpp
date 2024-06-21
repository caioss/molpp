#include "auxiliary.hpp"
#include "matchers.hpp"
#include <molpp/Atom.hpp>
#include <molpp/Residue.hpp>

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
    Topology& topology = data.topology();
    AtomData& atom_data = data.atoms();
    ResidueData& res_data = data.residues();
    res_data.resize(num_res);
    std::string const letters("ABCDEFGHIJKLMNOPQRSTUVWXYZ");

    // Set atoms
    for (index_t atom_idx = 0; atom_idx < num_atoms; atom_idx++)
    {
        std::string const code = letters.substr(atom_idx % 26, 1);

        index_t const res_idx = atom_idx / num_res_atoms;
        topology.link_entities({MolecularEntityCategory::Residue, res_idx}, {MolecularEntityCategory::Atom, atom_idx});
        atom_data.atomic_number(atom_idx) = atom_idx;
        atom_data.occupancy(atom_idx) = atom_idx;
        atom_data.temperature_factor(atom_idx) = atom_idx;
        atom_data.mass(atom_idx) = atom_idx;
        atom_data.charge(atom_idx) = atom_idx;
        atom_data.radius(atom_idx) = atom_idx;
        atom_data.name(atom_idx) = code;
        atom_data.type(atom_idx) = code;
        atom_data.alternate_location(atom_idx) = code;
        atom_data.insertion_code(atom_idx) = code;
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
        data.trajectory().add_timestep(Timestep(num_atoms));
        for (index_t atom_idx = 0; atom_idx < num_atoms; atom_idx++)
        {
            auto& coords = data.trajectory().timestep(frame_idx).coords();
            coords(Eigen::all, atom_idx) << atom_idx, atom_idx, atom_idx;
        }
    }

    return data;
}

TEST(Auxiliary, create_moldata) {
    MolData data = create_moldata(3, 2, 2, 1, 2);
    ASSERT_EQ(data.size<Atom>(), 6);

    // Atoms
    ASSERT_EQ(data.size<Atom>(), 6);
    std::vector<mol::Atom> atoms;
    for (index_t i = 0; i < data.size<Atom>(); ++i)
    {
        atoms.push_back(Atom(i, std::nullopt, data));
    }

    EXPECT_THAT(atoms, Pointwise(Prop(&Atom::residue_index),
                                 {0, 0, 1, 1, 2, 2}));
    EXPECT_THAT(atoms, Pointwise(Prop(&Atom::atomic_number),
                                 {0, 1, 2, 3, 4, 5}));
    EXPECT_THAT(atoms, Pointwise(PropFloat(&Atom::occupancy, 1e-5),
                                 {0, 1, 2, 3, 4, 5}));
    EXPECT_THAT(atoms, Pointwise(PropFloat(&Atom::temperature_factor, 1e-5),
                                 {0, 1, 2, 3, 4, 5}));
    EXPECT_THAT(atoms, Pointwise(PropFloat(&Atom::mass, 1e-5),
                                 {0, 1, 2, 3, 4, 5}));
    EXPECT_THAT(atoms, Pointwise(PropFloat(&Atom::charge, 1e-5),
                                 {0, 1, 2, 3, 4, 5}));
    EXPECT_THAT(atoms, Pointwise(PropFloat(&Atom::radius, 1e-5),
                                 {0, 1, 2, 3, 4, 5}));
    EXPECT_THAT(atoms, Pointwise(Prop(&Atom::name),
                                 {"A", "B", "C", "D", "E", "F"}));
    EXPECT_THAT(atoms, Pointwise(Prop(&Atom::type),
                                 {"A", "B", "C", "D", "E", "F"}));
    EXPECT_THAT(atoms, Pointwise(Prop(&Atom::alternate_location),
                                 {"A", "B", "C", "D", "E", "F"}));
    EXPECT_THAT(atoms, Pointwise(Prop(&Atom::insertion_code),
                                 {"A", "B", "C", "D", "E", "F"}));

    // Residues
    ASSERT_EQ(data.size<Residue>(), 3);
    std::vector<mol::Residue> residues;
    for (index_t i = 0; i < data.size<Residue>(); ++i)
    {
        residues.push_back(mol::Residue(i, std::nullopt, data));
    }
    EXPECT_THAT(residues, Pointwise(Prop(&Residue::id),
                                 {0, 1, 2}));
    EXPECT_THAT(residues, Pointwise(Prop(&Residue::name),
                                 {"A", "B", "C"}));

    // Bonds
    BondData const& bond_data = data.bonds();
    EXPECT_THAT(bond_data.bonded(0), UnorderedElementsAre(0, 2));
    EXPECT_THAT(bond_data.bonded(1), UnorderedElementsAre());
    EXPECT_THAT(bond_data.bonded(2), UnorderedElementsAre(0, 2, 4));
    EXPECT_THAT(bond_data.bonded(3), UnorderedElementsAre());
    EXPECT_THAT(bond_data.bonded(4), UnorderedElementsAre(2, 4));
    EXPECT_THAT(bond_data.bonded(5), UnorderedElementsAre());

    // Trajectory
    Trajectory const& traj_data = data.trajectory();
    EXPECT_EQ(traj_data.num_frames(), 2);
    EXPECT_THAT(traj_data.timestep(0).coords().reshaped(), ElementsAre(0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5));
    EXPECT_THAT(traj_data.timestep(1).coords().reshaped(), ElementsAre(0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5));
}
