#include <molpp/internal/SelIndex.hpp>
#include <numeric>
#include <algorithm>

using namespace mol::internal;

SelIndex::SelIndex(size_t const max_size)
: m_indices(max_size)
{
    std::iota(m_indices.begin(), m_indices.end(), 0);
}

bool SelIndex::contains(index_t const index) const
{
    auto const it = std::lower_bound(m_indices.begin(), m_indices.end(), index);
    return it != m_indices.end() && *it == index;
}
