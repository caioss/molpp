#include <molpp/internal/PropertyContainer.hpp>

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
}

mol::Property* PropertyContainer::emplace(property_key_type const key, bool const is_time_based, size_t const size, PropertyTrajectory::make_property_fn const make_property)
{
    auto [iter, result] = m_properties.emplace(key, PropertyTrajectory(is_time_based, make_property));

    if (!result)
    {
        return iter->second.get(0);
    }

    PropertyTrajectory& trajectory = iter->second;
    trajectory.add_frames(m_num_frames);
    trajectory.resize(size);
    return trajectory.get(0);
}

mol::Property const* PropertyContainer::get(property_key_type const key, Frame const frame) const
{
    auto iter = m_properties.find(key);
    if (iter == m_properties.end())
    {
        return nullptr;
    }
    return iter->second.get(frame.value_or(0));
}

mol::Property* PropertyContainer::get(property_key_type const key, Frame const frame)
{
    // Delegate to the const version
    return const_cast<mol::Property*>(std::as_const(*this).get(key, frame));
}
