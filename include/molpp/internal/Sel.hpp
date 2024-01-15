#ifndef SEL_HPP
#define SEL_HPP

#include <molpp/Property.hpp>
#include <molpp/MolError.hpp>
#include <molpp/MolppCore.hpp>
#include <molpp/internal/SelIndex.hpp>
#include <molpp/internal/requirements.hpp>
#include <molpp/internal/MolData.hpp>

#include <vector>
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
    {t.data_size(data)} -> std::same_as<size_t>;
    {t.atom_indices()} -> std::convertible_to<std::vector<index_t>>;
};

template <class Derived, class Other>
concept SelFromAtoms = requires(Derived sel, Other other, MolData data)
{
    {other.data()} -> std::same_as<MolData*>;
    {other.frame()} -> std::same_as<Frame>;
    {sel.from_atom_indices(other.atom_indices(), data)} -> std::same_as<SelIndex>;
};

template <class Type, class Derived>
class Sel
{
private:
    template <class ItType>
    class Iterator;

public:
    using value_type = Type;
    using iterator = Iterator<Type>;
    using const_iterator = Iterator<const Type>;
    using coords_type = Eigen::IndexedView<Coord3, Eigen::internal::AllRange<3>, std::vector<index_t>>;

    Sel() = delete;
    Sel(Sel &&) = default;
    Sel(Sel const&) = default;
    Sel& operator=(Sel &&) = default;
    Sel& operator=(Sel const&) = default;

    explicit Sel(SelIndex&& sel_index, MolData* data)
    requires SelDerived<Derived>
    : m_data{data}
    , m_index{sel_index}
    {
        if (m_data->properties().num_frames())
        {
            m_frame = 0;
        }
    }

    template <class Other>
    explicit Sel(Other&& other)
    requires SelDerived<Derived> && SelFromAtoms<Derived, Other>
    : Sel(Derived::from_atom_indices(other.atom_indices(), *(other.data())), other.data())
    {
        set_frame(other.frame());
    }

    explicit Sel(IndexRange auto const& indices, MolData* data)
    requires SelDerived<Derived>
    : Sel(SelIndex(indices, Derived::data_size(*data)), data)
    {}

    explicit Sel(MolData* data)
    requires SelDerived<Derived>
    : Sel(SelIndex(Derived::data_size(*data)), data)
    {}

    Frame frame() const
    {
        return m_frame;
    }

    void set_frame(Frame const frame)
    {
        if (frame && frame >= m_data->properties().num_frames())
        {
            throw mol::MolError("Out of bounds frame: " + std::to_string(*frame));
        }
        m_frame = frame;
    }

    size_t size() const
    {
        return m_index.size();
    }

    bool contains(index_t const index) const
    {
        return m_index.contains(index);
    }

    std::vector<index_t> const &indices() const
    {
        return m_index.indices();
    }

    iterator begin()
    {
        return iterator(m_data, m_index.indices_begin(), frame());
    }

    iterator end()
    {
        return iterator(m_data, m_index.indices_end(), frame());
    }

    const_iterator begin() const
    {
        return const_iterator(m_data, m_index.indices_begin(), frame());
    }

    const_iterator end() const
    {
        return const_iterator(m_data, m_index.indices_end(), frame());
    }

    Type operator[](size_t const index)
    {
        return Type(indices()[index], frame(), m_data);
    }

    Type at(size_t const index)
    {
        if (index >= size())
        {
            throw mol::MolError("Out of bounds index: " + std::to_string(index));
        }

        return Type(indices()[index], frame(), m_data);
    }

    Type by_index(size_t const index)
    {
        if (!contains(index))
        {
            throw mol::MolError("Atom index " + std::to_string(index) + " not found in the selection");
        }
        return Type(index, frame(), m_data);
    }

    Derived bonded()
    {
        Derived &derived = static_cast<Derived &>(*this);
        auto atom_indices = derived.atom_indices();
        std::vector<index_t> bonded_atoms = m_data->bonds().bonded(atom_indices.begin(), atom_indices.end());
        Derived sel(Derived::from_atom_indices(bonded_atoms, *m_data), m_data);
        sel.set_frame(frame());
        return sel;
    }

    std::vector<std::shared_ptr<mol::Bond>> bonds()
    {
        Derived &derived = static_cast<Derived &>(*this);
        auto atom_indices = derived.atom_indices();
        return m_data->bonds().bonds(atom_indices.begin(), atom_indices.end());
    }

    template <IsProperty PropertyType>
    bool has()
    {
        return property<PropertyType>();
    }

    template <IsProperty PropertyType>
    bool has() const
    {
        return property<PropertyType>();
    }

    template <IsProperty PropertyType>
    PropertyType* property()
    {
        return m_data->properties().template get<Type, PropertyType>(frame());
    }

    template <IsProperty PropertyType>
    PropertyType const* property() const
    {
        return m_data->properties().template get<Type, PropertyType>(frame());
    }

protected:
    MolData* data()
    {
        return m_data;
    };

    MolData const* data() const
    {
        return m_data;
    };


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

    Frame m_frame;
    MolData* m_data;
    SelIndex m_index;
};

} // namespace internal
} // namespace mol

#endif // SEL_HPP
