#ifndef SELINDEX_HPP
#define SELINDEX_HPP

#include <vector>
#include <numeric>

namespace mol::internal {

class SelIndex
{
public:
    using iterator = std::vector<size_t>::const_iterator;

    SelIndex() = delete;

    explicit SelIndex(size_t const max_size);
    explicit SelIndex(std::vector<size_t> const &indices, size_t const max_size);

    std::vector<size_t> const &indices() const
    {
        return m_indices;
    }

    std::vector<bool> const &selected() const
    {
        return m_selected;
    }

    size_t size() const
    {
        return m_indices.size();
    }

    iterator indices_begin() const
    {
        return m_indices.cbegin();
    }

    iterator indices_end() const
    {
        return m_indices.cend();
    }

    bool contains(size_t const index) const
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

private:
    std::vector<bool> m_selected;
    std::vector<size_t> m_indices;
};

} // namespace mol::internal

#endif // SELINDEX_HPP
