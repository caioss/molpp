#ifndef TIMESTEP_HPP
#define TIMESTEP_HPP

#include <Eigen/Dense>

namespace mol::internal
{

class Timestep
{
public:
    Timestep();
    Timestep(size_t const num_atoms);
    Timestep(const Timestep &src) = delete;
    Timestep &operator=(const Timestep &rhs) = delete;
    Timestep(Timestep &&src) noexcept;
    Timestep &operator=(Timestep &&rhs);
    void swap(Timestep &rhs);

    Eigen::Matrix3Xf &coords() { return m_coords; }

private:
    size_t m_num_atoms;
    Eigen::Matrix3Xf m_coords;
};

} // namespace mol::internal

#endif // TIMESTEP_HPP
