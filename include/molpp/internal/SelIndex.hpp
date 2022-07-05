#ifndef SELINDEX_HPP
#define SELINDEX_HPP

#include <molpp/MolError.hpp>
#include <molpp/MolppCore.hpp>
#include <vector>

namespace mol::internal {

class SelIndex
{
public:
    using value_type = size_t;
    using type = std::vector<value_type>;
    using iterator = std::vector<size_t>::const_iterator;

    SelIndex() = delete;

    explicit SelIndex(size_t const max_size);

    template <SelIndexCompatible T>
    explicit SelIndex(T const &indices, size_t const max_size)
    {
        // Find unique, sorted and in-bounds indices
        std::vector<bool> m_selected(max_size, false);
        size_t count = 0;
        for (size_t const& index : indices)
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

        // Fill m_indices with correct values
        m_indices.reserve(count);
        for (size_t i = 0; i < max_size; i++)
        {
            if (m_selected[i])
            {
                m_indices.push_back(i);
            }
        }
    }

    std::vector<size_t> const &indices() const
    {
        return m_indices;
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

    bool contains(size_t const index) const;

private:
    type m_indices;
};

} // namespace mol::internal

#endif // SELINDEX_HPP
