#include "AtomData.hpp"
#include "Atom.hpp"
#include <utility>

using namespace mol::internal;

AtomData::AtomData(size_t const num_atoms)
: m_index(num_atoms),
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
    for (size_t i = 0; i < num_atoms; ++i)
    {
        m_index[i] = i;
    }

}

std::shared_ptr<AtomData> AtomData::create(size_t const num_atoms)
{
    // Use *new* to access the private constructor
    return std::shared_ptr<AtomData>(new AtomData(num_atoms));
}

mol::Atom AtomData::index(size_t const index)
{
    return mol::Atom(index, shared_from_this());
}

void AtomData::add_timestep(Timestep &&ts)
{
    m_timestep.push_back(std::forward<Timestep &&>(ts));
}
