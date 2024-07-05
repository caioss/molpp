#include <molpp/internal/SelIndices.hpp>

#include <numeric>
#include <algorithm>

namespace mol::internal
{

SelIndices::SelIndices(size_t const max_size)
: m_sorted(true)
, m_indices(max_size)
{
    std::iota(m_indices.begin(), m_indices.end(), 0);
}

SelIndices::indices_type const& SelIndices::indices() const
{
    return m_indices;
}

SelIndices::value_type SelIndices::size() const
{
    return m_indices.size();
}

SelIndices::const_iterator SelIndices::begin() const
{
    return m_indices.cbegin();
}

SelIndices::const_iterator SelIndices::end() const
{
    return m_indices.cend();
}

bool SelIndices::contains(index_t const index) const
{
    if (m_sorted)
    {
        auto const it = std::lower_bound(m_indices.begin(), m_indices.end(), index);
        return it != m_indices.end() && *it == index;
    }
    else
    {
        return std::find(m_indices.begin(), m_indices.end(), index) != m_indices.end();
    }
}

} // namespace mol::internal
