#include <molpp/internal/SelIndex.hpp>
#include <molpp/MolError.hpp>

using namespace mol::internal;

SelIndex::SelIndex(size_t const max_size)
: m_selected(max_size, true),
  m_indices(max_size)
{
    std::iota(m_indices.begin(), m_indices.end(), 0);
}

SelIndex::SelIndex(std::vector<size_t> const &indices, size_t const max_size)
: m_selected(max_size, false)
{
    // Fill m_selected
    size_t count = 0;
    for (size_t index : indices)
    {
        if (index >= max_size)
        {
            throw mol::MolError("Out of bounds selection index: " + std::to_string(index));
        }

        if (!m_selected[index])
        {
            m_selected[index] = true;
            count++;
        }
    }

    // Fill m_indices with sorted and non-duplicated values
    m_indices.reserve(count);
    for (size_t i = 0; i < max_size; i++)
    {
        if (m_selected[i])
        {
            m_indices.push_back(i);
        }
    }
}
