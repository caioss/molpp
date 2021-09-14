#include "MolData.hpp"
#include <molpp/Atom.hpp>
#include <utility>

using namespace mol::internal;

MolData::MolData(size_t const num_atoms)
: m_num_atoms { num_atoms },
  m_properties(num_atoms),
  m_bonds(num_atoms)
{
}

std::shared_ptr<MolData> MolData::create(size_t const num_atoms)
{
    // Use *new* to access the private constructor
    return std::shared_ptr<MolData>(new MolData(num_atoms));
}
