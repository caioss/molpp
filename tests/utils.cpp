#include "utils.hpp"

using namespace mol;
using namespace mol::internal;

MolData create_moldata(size_t const num_res, size_t const num_res_atoms, size_t const num_chains, size_t const num_segments, size_t const num_frames)
{
    size_t const num_atoms { num_res * num_res_atoms };
    MolData data(num_atoms);
    Topology& topology = data.topology();
    AtomData& atom_data = data.atoms();
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
    ResidueData& res_data = data.residues();
    res_data.resize(num_res);
    if (num_res > 0)
    {
        topology.link_categories(MolecularEntityCategory::Residue, MolecularEntityCategory::Atom);
    }
    for (index_t res_idx = 0; res_idx < num_res; res_idx++)
    {
        res_data.id(res_idx) = res_idx;
        res_data.name(res_idx) = letters.substr(res_idx % 26, 1);;

        topology.link_entities({MolecularEntityCategory::Residue, res_idx}, {MolecularEntityCategory::Chain, res_idx % num_chains});

        topology.link_entities({MolecularEntityCategory::Residue, res_idx}, {MolecularEntityCategory::Segment, res_idx % num_segments});
    }

    // Set chains
    ChainData& chain_data = data.chains();
    chain_data.resize(num_chains);
    if (num_chains > 1)
    {
        topology.link_categories(MolecularEntityCategory::Chain, MolecularEntityCategory::Residue);
    }
    for (index_t chain_idx = 0; chain_idx < num_chains; chain_idx++)
    {
        chain_data.name(chain_idx) = letters.substr(chain_idx % 26, 1);
    }

    // Set segments
    SegmentData& segment_data = data.segments();
    segment_data.resize(num_segments);
    if (num_segments > 1)
    {
        topology.link_categories(MolecularEntityCategory::Segment, MolecularEntityCategory::Residue);
    }
    for (index_t segment_idx = 0; segment_idx < num_segments; segment_idx++)
    {
        segment_data.name(segment_idx) = letters.substr(segment_idx % 26, 1);
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
