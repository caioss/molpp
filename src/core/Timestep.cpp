#include "Timestep.hpp"

using namespace mol::internal;

Timestep::Timestep()
: m_num_atoms { 0 }
{}

Timestep::Timestep(size_t const num_atoms)
: m_num_atoms { num_atoms },
  m_coords(3, num_atoms)
{}

Timestep::Timestep(Timestep &&src) noexcept
: Timestep()
{
    swap(src);
}

Timestep &Timestep::operator=(Timestep &&rhs)
{
    if (this != &rhs)
    {
        swap(rhs);
    }
    return *this;
}

void Timestep::swap(Timestep &rhs)
{
    std::swap(this->m_num_atoms, rhs.m_num_atoms);
    std::swap(this->m_coords, rhs.m_coords);
}
