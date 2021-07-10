#include "AtomData.hpp"
#include "Atom.hpp"
#include <utility>

using namespace mol::internal;

AtomData::AtomData(size_t const num_atoms)
: m_num_atoms { num_atoms },
  m_properties(num_atoms),
  m_bonds(num_atoms)
{
}

std::shared_ptr<AtomData> AtomData::create(size_t const num_atoms)
{
    // Use *new* to access the private constructor
    return std::shared_ptr<AtomData>(new AtomData(num_atoms));
}

void AtomData::add_timestep(Timestep &&ts)
{
    m_timestep.push_back(std::forward<Timestep &&>(ts));
}
