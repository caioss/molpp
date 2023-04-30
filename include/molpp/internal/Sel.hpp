#ifndef SEL_HPP
#define SEL_HPP

#include <molpp/MolError.hpp>
#include <molpp/MolppCore.hpp>
#include <molpp/internal/BaseSel.hpp>
#include <molpp/internal/AtomAggregate.hpp>
#include <ranges>
#include <optional>
#include <concepts>

namespace mol {

class Bond;

namespace internal {

class MolData;

template <class Type, class Derived>
class Sel;

template <class Derived>
concept SelDerived = requires(Derived t, MolData data)
{
    std::derived_from<Derived, Sel<typename Derived::value_type, Derived>>;
    {t.atom_indices()} -> std::same_as<std::vector<index_t>>;
    {t.from_atom_indices({}, data)} -> std::same_as<std::vector<index_t>>;
    {t.max_size(data)} -> std::same_as<size_t>;
};

template <class Type>
concept SelConvertible = requires
{
    SelDerived<Type> || AtomAggregateDerived<Type>;
};

template <class Type, class Derived>
class Sel : public BaseSel
{
private:
    template <class ItType>
    class Iterator;

public:
    using value_type = Type;
    using iterator = Iterator<Type>;
    using const_iterator = Iterator<const Type>;
    using BaseSel::coords_type;

    Sel() = delete;

    template <SelConvertible Other>
    explicit Sel(Other&& other)
    requires SelDerived<Derived>
    : Sel(Derived::from_atom_indices(other.atom_indices(), *(other.data())), other.data())
    {
        set_frame(other.frame());
    }

    explicit Sel(MolData* data)
    requires SelDerived<Derived>
    : BaseSel(SelIndex(Derived::max_size(*data)), data)
    {}

    template <class T>
    explicit Sel(T const &indices, MolData* data)
    requires SelDerived<Derived> && SelIndexCompatible<T>
    : BaseSel(SelIndex(indices, Derived::max_size(*data)), data)
    {}

    iterator begin()
    {
        return iterator(data(), indices_begin(), frame());
    }

    iterator end()
    {
        return iterator(data(), indices_end(), frame());
    }

    const_iterator begin() const
    {
        return const_iterator(data(), indices_begin(), frame());
    }

    const_iterator end() const
    {
        return const_iterator(data(), indices_end(), frame());
    }

    Type operator[](size_t const index)
    {
        return Type(indices()[index], frame(), data());
    }

    Type at(size_t const index)
    {
        if (index >= size())
        {
            throw mol::MolError("Out of bounds index: " + std::to_string(index));
        }

        return Type(indices()[index], frame(), data());
    }

    Type by_index(size_t const index)
    {
        if (!contains(index))
        {
            throw mol::MolError("Atom index " + std::to_string(index) + " not found in the selection");
        }
        return Type(index, frame(), data());
    }

    coords_type coords()
    {
        Derived &derived = static_cast<Derived &>(*this);
        return coords(derived.atom_indices());
    }

    Derived bonded()
    {
        Derived &derived = static_cast<Derived &>(*this);
        Derived sel(Derived::from_atom_indices(bonded(derived.atom_indices()), *data()), data());
        sel.set_frame(frame());
        return sel;
    }

    std::vector<std::shared_ptr<mol::Bond>> bonds()
    {
        Derived &derived = static_cast<Derived &>(*this);
        return bonds(derived.atom_indices());
    }

protected:
    using BaseSel::coords;
    using BaseSel::bonded;
    using BaseSel::bonds;
    using BaseSel::data;

private:
    template <class ItType>
    class Iterator
    {
    private:
        using indices_iterator = SelIndex::iterator;

    public:
        using iterator_category = indices_iterator::iterator_category;
        using difference_type = indices_iterator::difference_type;
        using value_type = ItType;
        using pointer = ItType *;
        using reference = ItType &;

        Iterator(MolData* data, indices_iterator begin, Frame frame)
        : m_frame(frame),
          m_data(data),
          m_current(begin)
        {}

        value_type operator*() const
        {
            return ItType(*m_current, m_frame, m_data);
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
        Frame m_frame;
        MolData* m_data;
        indices_iterator m_current;
    };
};

} // namespace internal
} // namespace mol

#endif // SEL_HPP
