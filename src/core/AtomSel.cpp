#include "AtomSel.hpp"
#include "MolError.hpp"
#include "core/AtomData.hpp"
#include <numeric>
#include <algorithm>

using namespace mol;

AtomSel::AtomSel(std::shared_ptr<internal::AtomData> data)
: m_frame { 0 },
  m_selected(data->size(), true),
  m_indices(data->size()),
  m_data(data)
{
    std::iota(m_indices.begin(), m_indices.end(), 0);
}

AtomSel::AtomSel(std::vector<size_t> const &indices, std::shared_ptr<internal::AtomData> data)
: m_frame { 0 },
  m_selected(data->size(), false),
  m_indices(indices),
  m_data(data)
{
    update_indices();
}

AtomSel::AtomSel(std::vector<size_t> &&indices, std::shared_ptr<internal::AtomData> data)
: m_frame { 0 },
  m_selected(data->size(), false),
  m_indices(0),
  m_data(data)
{
    std::swap(m_indices, indices);
    update_indices();
}

void AtomSel::update_indices()
{
    std::sort(m_indices.begin(), m_indices.end());
    size_t const total_atoms = m_data->size();
    for (size_t index : m_indices)
    {
        if (index >= total_atoms)
        {
            throw mol::MolError("Out of bounds selection index: " + std::to_string(index));
        }

        m_selected[index] = true;
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
    return mol::Atom(m_indices[index], m_frame, m_data);
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
    if (index >= m_selected.size())
    {
        return false;
    }
    else
    {
        return m_selected[index];
    }
}

AtomSel::coords_type AtomSel::coords()
{
    return m_data->timestep(m_frame).coords()(Eigen::all, m_indices);
}

template<>
Atom AtomSel::Iterator<Atom>::operator*() const
{
    return mol::Atom(*m_current, m_frame, m_data);
}

template<>
const Atom AtomSel::Iterator<const Atom>::operator*() const
{
    return mol::Atom(*m_current, m_frame, m_data);
}
