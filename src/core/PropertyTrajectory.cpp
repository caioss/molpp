#include "PropertyTrajectory.hpp"

using namespace mol::internal;

PropertyTrajectory::PropertyTrajectory(bool const is_time_based, make_property_fn const make_property)
: m_make_property{make_property}
, m_time_based{is_time_based}
{
    if (!m_time_based)
    {
        emplace_frame();
    }
}

size_t PropertyTrajectory::is_time_based() const
{
    return m_time_based;
}

void PropertyTrajectory::resize(size_t const size)
{
    for (auto& property : m_trajectory)
    {
        property->resize(size);
    }
}

size_t PropertyTrajectory::num_frames() const
{
    return m_time_based ? m_trajectory.size() : 0;
}

mol::Property* PropertyTrajectory::add_frame()
{
    if (!m_time_based)
    {
        return nullptr;
    }

    return emplace_frame();
}

void PropertyTrajectory::remove_frame(Frame const frame)
{
    if (!m_time_based || !frame)
    {
        return;
    }

    auto iter = std::next(m_trajectory.begin(), frame.value());

    if (iter == m_trajectory.end())
    {
        return;
    }

    m_trajectory.erase(iter);
}

mol::Property* PropertyTrajectory::get(Frame const frame)
{
    auto iter = m_trajectory.begin();

    if (m_time_based && frame)
    {
        iter = std::next(iter, frame.value());
    }

    if (iter == m_trajectory.end())
    {
        return nullptr;
    }

    return iter->get();
}

mol::Property* PropertyTrajectory::emplace_frame()
{
    std::unique_ptr<Property>& property = m_trajectory.emplace_back();
    property = m_make_property();
    return property.get();
}
