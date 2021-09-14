#ifndef SEL_HPP
#define SEL_HPP

#include <molpp/internal/BaseSel.hpp>
#include <molpp/MolError.hpp>
#include <memory>
#include <optional>
#include <concepts>

namespace mol {

class Bond;
class AtomSel;
class ResidueSel;

namespace internal {

class MolData;

template <class Type, class Derived>
class Sel;

template <class Type, class Derived>
concept is_Sel = requires(Derived t)
{
    std::derived_from<Derived, Sel<Type, Derived>>;
    {t.atom_indices()} -> std::same_as<std::vector<size_t>>;
    {Derived::from_atom_indices({}, nullptr)} -> std::same_as<Derived>;
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

    Sel(std::shared_ptr<MolData> data)
    requires is_Sel<Type, Derived>
    : BaseSel(Derived::max_size(data), data)
    {}

    Sel(std::vector<size_t> const &indices, std::shared_ptr<MolData> data)
    requires is_Sel<Type, Derived>
    : BaseSel(Derived::max_size(data), indices, data)
    {}

    Sel(std::vector<size_t> &&indices, std::shared_ptr<MolData> data)
    requires is_Sel<Type, Derived>
    : BaseSel(Derived::max_size(data), std::forward<std::vector<size_t>&&>(indices), data)
    {}

    iterator begin()
    {
        return iterator(data(), indices_begin(), frame());
    }

    iterator end()
    {
        return iterator(data(), indices_end(), frame());
    }

    const_iterator cbegin()
    {
        return const_iterator(data(), indices_begin(), frame());
    }

    const_iterator cend()
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

    coords_type coords()
    {
        Derived &derived = static_cast<Derived &>(*this);
        return coords(derived.atom_indices());
    }

    std::shared_ptr<Derived> bonded()
    {
        Derived &derived = static_cast<Derived &>(*this);
        auto sel = std::make_shared<Derived>(Derived::from_atom_indices(bonded(derived.atom_indices()), data()));
        sel->set_frame(frame());
        return sel;
    }

    std::vector<std::shared_ptr<mol::Bond>> bonds()
    {
        Derived &derived = static_cast<Derived &>(*this);
        return bonds(derived.atom_indices());
    }

    std::shared_ptr<AtomSel> atoms()
    {
        Derived &derived = static_cast<Derived &>(*this);
        return atoms(derived.atom_indices());
    }

    std::shared_ptr<ResidueSel> residues()
    {
        Derived &derived = static_cast<Derived &>(*this);
        return residues(derived.atom_indices());
    }

protected:
    using BaseSel::coords;
    using BaseSel::bonded;
    using BaseSel::bonds;
    using BaseSel::atoms;
    using BaseSel::residues;

private:
    template <class ItType>
    class Iterator
    {
    private:
        using indices_iterator = std::vector<size_t>::iterator;

    public:
        using iterator_category = indices_iterator::iterator_category;
        using difference_type = indices_iterator::difference_type;
        using value_type = ItType;
        using pointer = ItType *;
        using reference = ItType &;

        Iterator(std::shared_ptr<MolData> data, indices_iterator begin, std::optional<size_t> frame)
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
        std::optional<size_t> m_frame;
        std::shared_ptr<MolData> m_data;
        indices_iterator m_current;
    };
};

} // namespace internal
} // namespace mol

#endif // SEL_HPP
