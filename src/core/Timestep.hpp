#ifndef TIMESTEP_HPP
#define TIMESTEP_HPP

#include <Eigen/Dense>

namespace mol::internal
{

class Timestep
{
public:
    Timestep()
    : m_num_atoms { 0 }
    {}

    Timestep(size_t const num_atoms)
    : m_num_atoms { num_atoms },
      m_coords(3, num_atoms)
    {}

    Timestep(const Timestep &src) = delete;
    Timestep &operator=(const Timestep &rhs) = delete;

    Timestep(Timestep &&src) noexcept
    : Timestep()
    {
        swap(src);
    }

    Timestep &operator=(Timestep &&rhs)
    {
        if (this != &rhs)
        {
            swap(rhs);
        }
        return *this;
    }

    void swap(Timestep &rhs)
    {
        std::swap(this->m_num_atoms, rhs.m_num_atoms);
        std::swap(this->m_coords, rhs.m_coords);
    }

    Eigen::Matrix3Xf &coords() { return m_coords; }

private:
    size_t m_num_atoms;
    Eigen::Matrix3Xf m_coords;
};

} // namespace mol::internal

#endif // TIMESTEP_HPP
