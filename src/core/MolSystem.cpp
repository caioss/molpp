#include "MolSystem.hpp"
#include "AtomData.hpp"
#include "readers/MolReader.hpp"
#include "MolError.hpp"

using namespace mol::internal;

mol::MolSystem::MolSystem(std::string const& topology)
{
    auto reader = MolReader::from_file(topology);
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
