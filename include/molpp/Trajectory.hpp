#ifndef TRAJECTORY_HPP
#define TRAJECTORY_HPP

#include <molpp/Timestep.hpp>
#include <vector>

namespace mol {

class Trajectory
{
public:
    Trajectory()
    {}

    size_t num_frames() const {
        return m_timestep.size();
    }

    Timestep &timestep(size_t const index) {
        return m_timestep[index];
    }

    Timestep const& timestep(size_t const index) const {
        return m_timestep[index];
    }

    void add_timestep(Timestep &&ts)
    {
        m_timestep.push_back(std::forward<Timestep>(ts));
    }

private:
    std::vector<Timestep> m_timestep;
};

} // namespace mol

#endif // TRAJECTORY_HPP
