#include <molpp/internal/PropertyTrajectory.hpp>

#include <utility>

using namespace mol::internal;

PropertyTrajectory::PropertyTrajectory(bool const is_time_based, make_property_fn const make_property)
: m_is_time_based{is_time_based}
, m_size{0}
, m_make_property{make_property}
{
    if (!m_is_time_based)
    {
        emplace_frame();
    }
}

size_t PropertyTrajectory::is_time_based() const
{
    return m_is_time_based;
}

void PropertyTrajectory::resize(size_t const size)
{
    if (size == m_size)
    {
        return;
    }

    for (auto& property : m_trajectory)
    {
        property->resize(size);
    }

    m_size = size;
}

size_t PropertyTrajectory::num_frames() const
{
    return m_is_time_based ? m_trajectory.size() : 0;
}

mol::Property* PropertyTrajectory::add_frame()
{
    if (!m_is_time_based)
    {
        return nullptr;
    }

    return emplace_frame();
}

void PropertyTrajectory::add_frames(size_t const num_frames)
{
    if (!m_is_time_based)
    {
        return;
    }

    for (size_t i = 0; i < num_frames; i++)
    {
        emplace_frame();
    }
}

void PropertyTrajectory::remove_frame(size_t const frame)
{
    if (!m_is_time_based)
    {
        return;
    }

    auto iter = std::next(m_trajectory.begin(), frame);

    if (iter == m_trajectory.end())
    {
        return;
    }

    m_trajectory.erase(iter);
}

mol::Property const* PropertyTrajectory::get(size_t const frame) const
{
    auto iter = m_trajectory.begin();

    if (m_is_time_based)
    {
        for (size_t i = 0; i < frame && iter != m_trajectory.end(); i++)
        {
            iter = std::next(iter);
        }
    }

    if (iter == m_trajectory.end())
    {
        return nullptr;
    }

    return iter->get();
}

mol::Property* PropertyTrajectory::get(size_t const frame)
{
    // Delegate to the const version
    return const_cast<mol::Property*>(std::as_const(*this).get(frame));
}

mol::Property* PropertyTrajectory::emplace_frame()
{
    std::unique_ptr<Property>& property = m_trajectory.emplace_back(m_make_property());
    property->resize(m_size);
    return property.get();
}
