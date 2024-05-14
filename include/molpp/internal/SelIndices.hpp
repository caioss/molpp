#ifndef MOLPP_INTERNAL_SELINDICES_HPP
#define MOLPP_INTERNAL_SELINDICES_HPP

#include <molpp/MolppCore.hpp>

#include <vector>
#include <unordered_set>

namespace mol::internal
{

class SelIndices
{
public:
    using value_type = index_t;
    using indices_type = std::vector<value_type>;
    using const_iterator = indices_type::const_iterator;

    SelIndices() = delete;
    explicit SelIndices(size_t const max_size);

    SelIndices(IndexRange auto const& indices)
    {
        // Remove duplicates and sort
        std::unordered_set<value_type> unique;
        for (value_type index : indices)
        {
            unique.insert(index);
        }
        m_indices.assign(unique.begin(), unique.end());
        std::sort(m_indices.begin(), m_indices.end());
    }

    indices_type const& indices() const;
    value_type size() const;
    const_iterator begin() const;
    const_iterator end() const;
    bool contains(index_t const index) const;

private:
    indices_type m_indices;
};

} // namespace mol::internal

#endif // MOLPP_INTERNAL_SELINDICES_HPP
