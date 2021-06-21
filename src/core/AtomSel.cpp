#include "AtomSel.hpp"
#include "MolError.hpp"
#include "core/AtomData.hpp"
#include <numeric>
#include <algorithm>

using namespace mol;

AtomSel::AtomSel(std::shared_ptr<internal::AtomData> data)
: m_frame { 0 },
  m_indices(data->size()),
  m_data(data)
{
    std::iota(m_indices.begin(), m_indices.end(), 0);
}

AtomSel::AtomSel(std::vector<size_t> const &indices, std::shared_ptr<internal::AtomData> data)
: m_frame { 0 },
  m_data(data)
{
    // Remove duplicates
    size_t num_atoms = 0;
    size_t const total_atoms = data->size();
    std::vector<bool> selected(total_atoms, false);
    for (size_t index : indices)
    {
        if (index >= total_atoms)
        {
            throw mol::MolError("Out of bounds selection index: " + std::to_string(index));
        }
        if (!selected[index])
        {
            ++num_atoms;
            selected[index] = true;
        }
    }

    // Fill with sorted indices
    m_indices.reserve(num_atoms);
    for (size_t index = 0; index < total_atoms; ++index)
    {
        if (selected[index])
        {
            m_indices.push_back(index);
        }
    }
}

AtomSel::iterator AtomSel::begin()
{
    return AtomSel::iterator(m_data, m_indices.begin(), m_frame);
}

AtomSel::iterator AtomSel::end()
{
    return AtomSel::iterator(m_data, m_indices.end(), m_frame);
}

AtomSel::const_iterator AtomSel::cbegin()
{
    return AtomSel::const_iterator(m_data, m_indices.begin(), m_frame);
}

AtomSel::const_iterator AtomSel::cend()
{
    return AtomSel::const_iterator(m_data, m_indices.end(), m_frame);
}

Atom AtomSel::operator[](size_t const index)
{
    return m_data->index(index, m_frame);
}

void AtomSel::set_frame(size_t frame)
{
    if (frame >= m_data->num_frames())
    {
        throw mol::MolError("Out of bounds selection frame: " + std::to_string(frame));
    }
    m_frame = frame;
}

bool AtomSel::contains(size_t const index) const
{
    return std::binary_search(m_indices.begin(), m_indices.end(), index);
}

AtomSel::coords_type AtomSel::coords()
{
    return m_data->timestep(m_frame).coords()(Eigen::all, m_indices);
}

template<>
Atom AtomSel::Iterator<Atom>::operator*() const
{
    return m_data->index(*m_current, m_frame);
}

template<>
const Atom AtomSel::Iterator<const Atom>::operator*() const
{
    return m_data->index(*m_current, m_frame);
}
