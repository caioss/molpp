#include "MolData.hpp"

using namespace mol::internal;

MolData::MolData(size_t const num_atoms)
: m_num_atoms{num_atoms},
  m_properties_old(num_atoms),
  m_bonds(num_atoms)
{
}
