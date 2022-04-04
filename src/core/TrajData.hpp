#ifndef TRAJDATA_HPP
#define TRAJDATA_HPP

#include "core/Timestep.hpp"
#include <vector>

namespace mol::internal {

class TrajData
{
public:
    TrajData()
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

} // namespace mol::internal

#endif // TRAJDATA_HPP
