#include "AtomData.hpp"
#include "Atom.hpp"
#include <utility>

using namespace mol::internal;

AtomData::AtomData(size_t const num_atoms)
: m_num_atoms { num_atoms },
  m_resid(num_atoms, -1),
  m_residue(num_atoms, -1),
  m_atomic(num_atoms, 0),
  m_occupancy(num_atoms, 0),
  m_tempfactor(num_atoms, 0),
  m_mass(num_atoms, 0),
  m_charge(num_atoms, 0),
  m_radius(num_atoms, 0),
  m_name(num_atoms),
  m_type(num_atoms),
  m_resname(num_atoms),
  m_segid(num_atoms),
  m_chain(num_atoms),
  m_altloc(num_atoms)
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
