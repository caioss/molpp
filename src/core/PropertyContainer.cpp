#include "core/PropertyContainer.hpp"

using namespace mol::internal;

PropertyContainer::PropertyContainer()
: m_num_frames{0}
{
}

size_t PropertyContainer::num_frames() const
{
    return m_num_frames;
}

mol::Frame PropertyContainer::add_frame()
{
    for (auto& [key, trajectory] : m_properties)
    {
        trajectory.add_frame();
    }

    return m_num_frames++;
}

void PropertyContainer::remove_frame(Frame const frame)
{
    if (!frame || frame >= m_num_frames)
    {
        return;
    }

    for (auto& [key, trajectory] : m_properties)
    {
        trajectory.remove_frame(frame.value());
    }

    m_num_frames--;
}

size_t PropertyContainer::size(size_key_type const key) const
{
    auto const iter = m_sizes.find(key);
    if (iter == m_sizes.end())
    {
        return 0;
    }
    return iter->second;
}

void PropertyContainer::resize(size_key_type const key, size_t const size)
{
    for (auto& [property_key, trajectory] : m_properties)
    {
        if (property_key.first != key)
        {
            continue;
        }

        trajectory.resize(size);
    }
    m_sizes[key] = size;
}

mol::Property* PropertyContainer::emplace(property_key_type const key, bool const is_time_based, PropertyTrajectory::make_property_fn const make_property)
{
    auto [iter, result] = m_properties.emplace(key, PropertyTrajectory(is_time_based, make_property));

    if (!result)
    {
        return iter->second.get(0);
    }

    PropertyTrajectory& trajectory = iter->second;
    trajectory.add_frames(m_num_frames);
    trajectory.resize(size(key.first));
    return trajectory.get(0);
}

mol::Property* PropertyContainer::get(property_key_type const key, Frame const frame)
{
    auto iter = m_properties.find(key);
    if (iter == m_properties.end())
    {
        return nullptr;
    }
    return iter->second.get(frame.value_or(0));
}
