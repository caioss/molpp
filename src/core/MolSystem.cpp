#include "MolSystem.hpp"
#include "AtomData.hpp"
#include "AtomSel.hpp"
#include "MolError.hpp"
#include "readers/MolReader.hpp"
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

    m_atoms = reader->read_topology(topology);
    if (!m_atoms)
    {
        throw mol::MolError("Error reading file " + topology);
    }
}
