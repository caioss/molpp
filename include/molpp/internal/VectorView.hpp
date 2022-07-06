#ifndef VECTORVIEW_HPP
#define VECTORVIEW_HPP

#include <ranges>
#include <iterator>

namespace mol {

template <std::ranges::random_access_range Container, std::ranges::bidirectional_range Index>
class SequenceView
{
private:
    template <class ItValueType, std::bidirectional_iterator ItType>
    class Iterator;

public:
    using value_type = typename Container::value_type;
    using pointer = value_type*;
    using reference = value_type&;
    using index_type = typename Index::value_type;
    using iterator = Iterator<value_type, typename Index::iterator>;
    using const_iterator = Iterator<value_type const, typename Index::const_iterator>;

    SequenceView(Container& data, Index const& indices)
    : m_data{data},
      m_indices{indices}
    {}

    size_t size() const
    {
        return m_indices.size();
    }

    iterator begin()
    {
        return iterator(m_data, m_indices.begin());
    }

    iterator end()
    {
        return iterator(m_data, m_indices.end());
    }

    const_iterator begin() const
    {
        return const_iterator(m_data, m_indices.cbegin());
    }

    const_iterator end() const
    {
        return const_iterator(m_data, m_indices.cend());
    }

    reference operator[](size_t const index)
    {
        return m_data[m_indices[index]];
    }

    reference at(size_t const index)
    {
        return m_data.at(m_indices.at(index));
    }

private:
    Container& m_data;
    Index const& m_indices;

    template <class ItValueType, std::bidirectional_iterator ItType>
    class Iterator
    {
    public:
        using iterator_category = ItType::iterator_category;
        using difference_type = ItType::difference_type;
        using value_type = ItValueType;
        using pointer = value_type*;
        using reference = value_type&;

        Iterator(Container& data, ItType begin)
        : m_data{data},
          m_current{begin}
        {}

        reference operator*() const
        {
            return m_data[*m_current];
        }

        Iterator& operator++()
        {
            m_current++;
            return *this;
        }

        Iterator operator++(int)
        {
            Iterator tmp = *this;
            ++(*this);
            return tmp;
        }

        difference_type operator-(const Iterator& other)
        {
            return m_current - other.m_current;
        }

        bool operator==(const Iterator& other)
        {
            return m_current == other.m_current;
        };

        bool operator!=(const Iterator& other)
        {
            return m_current != other.m_current;
        };

    private:
        Container& m_data;
        ItType m_current;
    };
};

} // namespace mol

#endif // VECTORVIEW_HPP
