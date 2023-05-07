#include <molpp/AtomSel.hpp>
#include <molpp/MolError.hpp>
#include <molpp/MolSystem.hpp>
#include <molpp/ResidueSel.hpp>
#include "core/MolData.hpp"
#include "readers/MolReader.hpp"
#include "guessers/AtomBondGuesser.hpp"
#include "guessers/ResidueBondGuesser.hpp"
#include <filesystem>

using namespace mol;
using namespace mol::internal;

MolSystem::MolSystem(std::string const &topology)
{
    auto reader = MolReader::from_file_ext(std::filesystem::path(topology).extension());
    if (!reader)
    {
        throw mol::MolError("No reader for file " + topology);
    }

    m_data = reader->read_topology(topology);
    if (!m_data)
    {
        throw mol::MolError("Error reading file " + topology);
    }
}

mol::MolSystem::~MolSystem()
{
}

void MolSystem::add_trajectory(std::string const& file_name, int begin, int end, int step)
{
    auto reader = MolReader::from_file_ext(std::filesystem::path(file_name).extension());
    if (!reader)
    {
        throw mol::MolError("No reader for file " + file_name);
    }

    MolReader::Status status = reader->read_trajectory(file_name, *m_data, begin, end, step);

    if (status != MolReader::SUCCESS)
    {
        switch (status)
        {
            case MolReader::WRONG_ATOMS:
                throw mol::MolError("Trajectory with wrong number of atoms");

            default:
                throw mol::MolError("Error reading file " + file_name);
        }
    }
}

AtomSel MolSystem::atoms(Frame const frame) const
{
    AtomSel sel(m_data.get());
    sel.set_frame(frame);
    return sel;
}

AtomSel MolSystem::select(std::vector<index_t> const& indices, Frame const frame) const
{
    AtomSel sel(indices, m_data.get());
    sel.set_frame(frame);
    return sel;
}

AtomSel MolSystem::select(std::string const& selection, Frame const frame) const
{
    if (selection == "all")
    {
        return atoms(frame);
    }

    return selector(selection).apply(frame);
}

AtomSelector MolSystem::selector(std::string const& selection) const
{
    return AtomSelector(selection, m_data.get());
}

void MolSystem::reset_bonds()
{
    m_data->bonds().clear();
}

void MolSystem::guess_bonds(Frame const frame)
{
    AtomSel all_atoms = atoms(frame);
    ResidueSel all_residues(all_atoms);

    // Fill tabulated bonds first
    ResidueBondGuesser res_guesser;
    res_guesser.apply(all_residues);

    // Bond heuristics are the last step
    if (m_data->trajectory().num_frames())
    {
        AtomBondGuesser atom_guesser;
        atom_guesser.apply(all_atoms);
    }
}
