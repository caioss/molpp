#ifndef TIMESTEP_HPP
#define TIMESTEP_HPP

#include <Eigen/Dense>

namespace mol::internal
{

class Timestep
{
public:
    using coords_type = Eigen::Matrix3Xf;

    Timestep();
    Timestep(size_t const num_atoms);
    Timestep(const Timestep &src) = delete;
    Timestep &operator=(const Timestep &rhs) = delete;
    Timestep(Timestep &&src) noexcept;
    Timestep &operator=(Timestep &&rhs);
    void swap(Timestep &rhs);

    coords_type &coords() { return m_coords; }

private:
    size_t m_num_atoms;
    coords_type m_coords;
};

} // namespace mol::internal

#endif // TIMESTEP_HPP
