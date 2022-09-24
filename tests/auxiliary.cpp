#include "auxiliary.hpp"
#include "matchers.hpp"
#include <molpp/Atom.hpp>
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <optional>

using namespace testing;
using namespace mol;
using namespace mol::internal;

std::shared_ptr<MolData> create_moldata(size_t const num_res, size_t const num_res_atoms, size_t const num_chains, size_t const num_segments, size_t const num_frames)
{
    size_t const num_atoms { num_res * num_res_atoms };
    auto data = MolData::create(num_atoms);
    AtomData& atom_data = data->atoms();
    ResidueData& res_data = data->residues();
    data->residues().resize(num_res);
    std::string const letters("ABCDEFGHIJKLMNOPQRSTUVWXYZ");

    // Set atoms
    for (index_t atom_idx = 0; atom_idx < num_atoms; atom_idx++)
    {
        std::string const code = letters.substr(atom_idx % 26, 1);

        index_t const res_idx = atom_idx / num_res_atoms;
        atom_data.residue(atom_idx) = res_idx;
        res_data.add_atom(res_idx, atom_idx);
        atom_data.atomic(atom_idx) = atom_idx;
        atom_data.occupancy(atom_idx) = atom_idx;
        atom_data.tempfactor(atom_idx) = atom_idx;
        atom_data.mass(atom_idx) = atom_idx;
        atom_data.charge(atom_idx) = atom_idx;
        atom_data.radius(atom_idx) = atom_idx;
        atom_data.name(atom_idx) = code;
        atom_data.type(atom_idx) = code;
        atom_data.altloc(atom_idx) = code;
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
        data->bonds().add_bond(atom_idx, atom_idx + num_res_atoms);
    }

    // Trajectory
    for (size_t frame_idx = 0; frame_idx < num_frames; frame_idx++)
    {
        data->trajectory().add_timestep(Timestep(num_atoms));
        for (index_t atom_idx = 0; atom_idx < num_atoms; atom_idx++)
        {
            auto& coords = data->trajectory().timestep(frame_idx).coords();
            coords(Eigen::all, atom_idx) << atom_idx, atom_idx, atom_idx;
        }
    }

    return data;
}

TEST(Auxiliary, create_moldata) {
    auto data = create_moldata(3, 2, 2, 1, 2);
    ASSERT_THAT(data, NotNull());
    EXPECT_EQ(data->size(), 6);

    // Atoms
    EXPECT_EQ(data->atoms().size(), 6);
    std::vector<mol::Atom> atoms;
    for (index_t i = 0; i < data->size(); ++i)
    {
        atoms.push_back(Atom(i, std::nullopt, data));
    }

    EXPECT_THAT(atoms, Pointwise(Prop(&Atom::residue_id),
                                 {0, 0, 1, 1, 2, 2}));
    EXPECT_THAT(atoms, Pointwise(Prop(&Atom::atomic),
                                 {0, 1, 2, 3, 4, 5}));
    EXPECT_THAT(atoms, Pointwise(PropFloat(&Atom::occupancy, 1e-5),
                                 {0, 1, 2, 3, 4, 5}));
    EXPECT_THAT(atoms, Pointwise(PropFloat(&Atom::tempfactor, 1e-5),
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
    EXPECT_THAT(atoms, Pointwise(Prop(&Atom::altloc),
                                 {"A", "B", "C", "D", "E", "F"}));
    EXPECT_THAT(atoms, Pointwise(Prop(&Atom::resid),
                                 {0, 0, 1, 1, 2, 2}));
    EXPECT_THAT(atoms, Pointwise(Prop(&Atom::resname),
                                 {"A", "A", "B", "B", "C", "C"}));
    EXPECT_THAT(atoms, Pointwise(Prop(&Atom::segid),
                                 {"A", "A", "A", "A", "A", "A"}));
    EXPECT_THAT(atoms, Pointwise(Prop(&Atom::chain),
                                 {"A", "A", "B", "B", "A", "A"}));

    // Residues
    EXPECT_EQ(data->residues().size(), 3);

    // Bonds
    BondData const& bond_data = data->bonds();
    EXPECT_THAT(bond_data.bonded(0), UnorderedElementsAre(0, 2));
    EXPECT_THAT(bond_data.bonded(1), UnorderedElementsAre());
    EXPECT_THAT(bond_data.bonded(2), UnorderedElementsAre(0, 2, 4));
    EXPECT_THAT(bond_data.bonded(3), UnorderedElementsAre());
    EXPECT_THAT(bond_data.bonded(4), UnorderedElementsAre(2, 4));
    EXPECT_THAT(bond_data.bonded(5), UnorderedElementsAre());

    // Trajectory
    Trajectory const& traj_data = data->trajectory();
    EXPECT_EQ(traj_data.num_frames(), 2);
    EXPECT_THAT(traj_data.timestep(0).coords().reshaped(), ElementsAre(0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5));
    EXPECT_THAT(traj_data.timestep(1).coords().reshaped(), ElementsAre(0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5));
}
