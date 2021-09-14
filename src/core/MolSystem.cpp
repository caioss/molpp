#include "MolData.hpp"
#include "readers/MolReader.hpp"
#include <molpp/MolSystem.hpp>
#include <molpp/AtomSel.hpp>
#include <molpp/MolError.hpp>
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

void MolSystem::add_trajectory(std::string const &file_name, int begin, int end, int step)
{
    auto reader = MolReader::from_file_ext(std::filesystem::path(file_name).extension());
    if (!reader)
    {
        throw mol::MolError("No reader for file " + file_name);
    }

    MolReader::Status status = reader->read_trajectory(file_name, m_data, begin, end, step);

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

std::shared_ptr<AtomSel> MolSystem::atoms() const
{
    return std::make_shared<AtomSel>(m_data);
}

std::shared_ptr<AtomSel> MolSystem::select(std::vector<size_t> const &indices) const
{
    return std::make_shared<AtomSel>(indices, m_data);
}
