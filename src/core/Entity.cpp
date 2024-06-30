#include <molpp/internal/Entity.hpp>

namespace mol
{

internal::Entity::Entity(index_t const index, Frame const frame, internal::MolData& data)
: m_index{index}
, m_frame(frame)
, m_data{&data}
{}

bool internal::Entity::operator==(Entity const& other) const
{
    return m_data == other.m_data && m_index == other.m_index && m_frame == other.m_frame;
}

index_t internal::Entity::index() const
{
    return m_index;
}

std::vector<index_t> internal::Entity::indices() const
{
    return std::vector<index_t>{m_index};
}

Frame internal::Entity::frame() const
{
    return m_frame;
}

void internal::Entity::set_frame(Frame const frame)
{
    if (frame && frame >= m_data->trajectory().num_frames())
    {
        throw Error("Out of bounds frame: " + std::to_string(*frame));
    }
    m_frame = frame;
}

internal::MolData& internal::Entity::data()
{
    return *m_data;
}

internal::MolData const& internal::Entity::data() const
{
    return *m_data;
}

} // namespace mol
