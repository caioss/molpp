#ifndef MOLPP_TIMESTEP_HPP
#define MOLPP_TIMESTEP_HPP

#include <molpp/MolppCore.hpp>

namespace mol
{

class Timestep
{
public:
    Timestep();
    Timestep(size_t const num_atoms);
    Timestep(Timestep const& src) = delete;
    Timestep& operator=(Timestep const& rhs) = delete;
    Timestep(Timestep&& src) noexcept;
    Timestep& operator=(Timestep&& rhs);
    void swap(Timestep& rhs);

    Coord3& coords()
    {
        return m_coords;
    }

    Coord3 const& coords() const
    {
        return m_coords;
    }

private:
    size_t m_num_atoms;
    Coord3 m_coords;
};

} // namespace mol

#endif // MOLPP_TIMESTEP_HPP
