#ifndef ATOMSEL_HPP
#define ATOMSEL_HPP

#include "Atom.hpp"
#include "MolppCore.hpp"
#include <vector>
#include <memory>

namespace mol {

namespace internal {
class AtomData;
}

class AtomSel
{
template <class Type>
class Iterator;

public:
    using value_type = Atom;
    using iterator = Iterator<Atom>;
    using const_iterator = Iterator<const Atom>;
    using coords_type = Eigen::IndexedView<Coord3, Eigen::internal::AllRange<3>, std::vector<size_t>>;

    AtomSel() = delete;
    AtomSel(std::shared_ptr<internal::AtomData> data);
    AtomSel(std::vector<size_t> const &indices, std::shared_ptr<internal::AtomData> data);

    iterator begin();
    iterator end();
    const_iterator cbegin();
    const_iterator cend();
    Atom operator[](size_t const index);
    size_t size() const { return m_indices.size(); }

    size_t frame() const { return m_frame; }
    void set_frame(size_t frame);
    bool contains(size_t const index) const;
    std::vector<size_t> const &indices() const { return m_indices; }
    coords_type coords();

private:
    size_t m_frame;
    std::vector<bool> m_selected;
    std::vector<size_t> m_indices;
    std::shared_ptr<internal::AtomData> m_data;

template <class Type>
class Iterator
{
private:
    using indices_iterator = std::vector<size_t>::iterator;

public:
    using iterator_category = indices_iterator::iterator_category;
    using difference_type = indices_iterator::difference_type;
    using value_type = Type;
    using pointer = Type *;
    using reference = Type &;

    Iterator(std::shared_ptr<internal::AtomData> data, indices_iterator begin, size_t frame)
    : m_frame { frame },
      m_data(data),
      m_current(begin)
    {}
    value_type operator*() const;
    Iterator& operator++() { m_current++; return *this; }
    Iterator operator++(int) { Iterator tmp = *this; ++(*this); return tmp; }
    difference_type operator-(const Iterator& other) { return m_current - other.m_current; }
    bool operator==(const Iterator& other) { return m_current == other.m_current; };
    bool operator!=(const Iterator& other) { return m_current != other.m_current; };

private:
    size_t m_frame;
    std::shared_ptr<internal::AtomData> m_data;
    indices_iterator m_current;
};

};

} // namespace mol

#endif // ATOMSEL_HPP
