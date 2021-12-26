#ifndef ITERATORS_HPP
#define ITERATORS_HPP

namespace mol::internal {

template <class Iterator>
class IteratorWrapper
{
public:
    using iterator_category = typename Iterator::iterator_category;
    using difference_type = typename Iterator::difference_type;
    using value_type = typename Iterator::value_type;
    using pointer = value_type *;
    using reference = value_type &;

    IteratorWrapper(Iterator iterator)
    : m_iterator(iterator)
    {}

    IteratorWrapper& operator++()
    {
        m_iterator++;
        return *this;
    }

    IteratorWrapper operator++(int)
    {
        IteratorWrapper tmp = *this;
        ++(*this);
        return tmp;
    }

    difference_type operator-(const IteratorWrapper& other)
    {
        return m_iterator - other.m_iterator;
    }

    bool operator==(IteratorWrapper const &other) const
    {
        return m_iterator == other.m_iterator;
    };

    bool operator!=(IteratorWrapper const &other) const
    {
        return !(m_iterator == other.m_iterator);
    };

protected:
    Iterator m_iterator;
};

} // namespace mol::internal

#endif // ITERATORS_HPP
