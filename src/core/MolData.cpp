#include <molpp/internal/MolData.hpp>

using namespace mol::internal;

MolData::MolData(size_t const num_atoms)
: m_properties_old(num_atoms)
, m_bonds(num_atoms)
{
}

size_t MolData::num_frames() const
{
    return m_properties.num_frames();
}

mol::Frame MolData::add_frame()
{
    return m_properties.add_frame();
}

void MolData::remove_frame(Frame const frame)
{
    m_properties.remove_frame(frame);
}
