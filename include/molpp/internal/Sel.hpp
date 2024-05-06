#ifndef MOLPP_INTERNAL_SEL_HPP
#define MOLPP_INTERNAL_SEL_HPP

#include <molpp/MolError.hpp>
#include <molpp/MolppCore.hpp>
#include <molpp/internal/SelIndices.hpp>
#include <molpp/internal/MolData.hpp>

#include <vector>
#include <concepts>

namespace mol
{

class Bond;

namespace internal
{

class MolData;

template<class Type, class Derived>
class Sel;

template<class LHS, class RHS>
concept SelConvertible = requires(LHS lhs, RHS rhs) {
    { rhs.data() } -> std::same_as<MolData*>;
    { rhs.frame() } -> std::same_as<Frame>;
    { rhs.as_atom_indices() } -> internal::IndexRange;
    { lhs.from_atom_indices(rhs.as_atom_indices(), std::declval<MolData>()) } -> internal::IndexRange;
};

template<class Type, class Derived>
class Sel
{
private:
    template<class ItType>
    class Iterator;

public:
    using indices_type = std::vector<index_t>;
    using value_type = Type;
    using iterator = Iterator<Type>;
    using const_iterator = Iterator<const Type>;

    Sel() = delete;
    Sel(Sel&&) = default;
    Sel(Sel const&) = default;
    Sel& operator=(Sel&&) = default;
    Sel& operator=(Sel const&) = default;

    explicit Sel(SelIndices&& sel_index, MolData* data)
    : m_data{data}
    , m_index{std::forward<SelIndices>(sel_index)}
    {
        if (m_data->trajectory().num_frames())
        {
            m_frame = 0;
        }
    }

    template<class RHS>
    explicit Sel(RHS&& rhs)
    requires SelConvertible<Derived, RHS>
    : Sel(SelIndices(Derived::from_atom_indices(rhs.as_atom_indices(), *(rhs.data()))), rhs.data())
    {
        set_frame(rhs.frame());
    }

    explicit Sel(IndexRange auto const& indices, MolData* data)
    : Sel(SelIndices(indices), data)
    {}

    explicit Sel(MolData* data)
    : Sel(SelIndices(data->size<Type>()), data)
    {}

    Frame frame() const
    {
        return m_frame;
    }

    void set_frame(Frame const frame)
    {
        if (frame && frame >= m_data->trajectory().num_frames())
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

    std::vector<index_t> const& indices() const
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
    template<class ItType>
    class Iterator
    {
    private:
        using indices_iterator = SelIndices::const_iterator;

    public:
        using iterator_category = indices_iterator::iterator_category;
        using difference_type = indices_iterator::difference_type;
        using value_type = ItType;
        using pointer = ItType*;
        using reference = ItType&;

        Iterator(MolData* data, indices_iterator begin, Frame frame)
        : m_frame(frame)
        , m_data(data)
        , m_current(begin)
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
    SelIndices m_index;
};

} // namespace internal
} // namespace mol

#endif // MOLPP_INTERNAL_SEL_HPP
