#include <molpp/internal/MolecularEntity.hpp>

namespace mol
{

internal::MolecularEntity::MolecularEntity(index_t const index, Frame const frame, internal::MolData* data)
: m_index{index}
, m_frame(frame)
, m_data{data}
{}

bool internal::MolecularEntity::operator==(MolecularEntity const& other) const
{
    return m_data == other.m_data && m_index == other.m_index && m_frame == other.m_frame;
}

index_t internal::MolecularEntity::index() const
{
    return m_index;
}

Frame internal::MolecularEntity::frame() const
{
    return m_frame;
}

void internal::MolecularEntity::set_frame(Frame const frame)
{
    if (frame && frame >= m_data->trajectory().num_frames())
    {
        throw MolError("Out of bounds frame: " + std::to_string(*frame));
    }
    m_frame = frame;
}

bool internal::MolecularEntity::is_valid() const
{
    return m_data;
}

internal::MolData* internal::MolecularEntity::data()
{
    return m_data;
}

internal::MolData const* internal::MolecularEntity::data() const
{
    return m_data;
}

} // namespace mol
