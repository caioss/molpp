#ifndef MOLPP_INTERNAL_SEL_HPP
#define MOLPP_INTERNAL_SEL_HPP

#include <molpp/Error.hpp>
#include <molpp/Common.hpp>
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

template<class Derived>
class Sel;

template<class RHS>
concept SelConvertible = requires(RHS rhs) {
    { rhs.data() } -> std::same_as<MolData&>;
    { rhs.frame() } -> std::same_as<Frame>;
    { rhs.category() } -> std::same_as<EntityCategory>;
    { rhs.indices() } -> internal::IndexRange;
};

template<class Derived>
struct SelTraits;

template<class Derived>
class Sel
{
private:
    template<class ItType>
    class Iterator;

public:
    using indices_type = std::vector<index_t>;
    using value_type = typename SelTraits<Derived>::entity_type;
    using iterator = Iterator<value_type>;
    using const_iterator = Iterator<const value_type>;

    Sel() = delete;
    Sel(Sel&&) = default;
    Sel(Sel const&) = default;
    Sel& operator=(Sel&&) = default;
    Sel& operator=(Sel const&) = default;

    explicit Sel(SelIndices&& sel_index, MolData& data)
    : m_data{&data}
    , m_index{std::forward<SelIndices>(sel_index)}
    {
        if (m_data->trajectory().num_frames())
        {
            m_frame = 0;
        }
    }

    template<class RHS>
    explicit Sel(RHS&& rhs)
    requires SelConvertible<RHS>
    : Sel(SelIndices(rhs.data().topology().convert(rhs.indices(), rhs.category(), category())), rhs.data())
    {
        set_frame(rhs.frame());
    }

    explicit Sel(IndexRange auto const& indices, MolData& data)
    : Sel(SelIndices(indices), data)
    {}

    explicit Sel(MolData& data)
    : Sel(SelIndices(data.size<value_type>()), data)
    {}

    static EntityCategory category()
    {
        return value_type::category();
    }

    Frame frame() const
    {
        return m_frame;
    }

    void set_frame(Frame const frame)
    {
        if (frame && frame >= m_data->trajectory().num_frames())
        {
            throw mol::Error("Out of bounds frame: " + std::to_string(*frame));
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
        return iterator(*m_data, m_index.begin(), frame());
    }

    iterator end()
    {
        return iterator(*m_data, m_index.end(), frame());
    }

    const_iterator begin() const
    {
        return const_iterator(*m_data, m_index.begin(), frame());
    }

    const_iterator end() const
    {
        return const_iterator(*m_data, m_index.end(), frame());
    }

    value_type operator[](size_t const index)
    {
        return value_type(indices()[index], frame(), *m_data);
    }

    value_type at(size_t const index)
    {
        if (index >= size())
        {
            throw mol::Error("Out of bounds index: " + std::to_string(index));
        }

        return value_type(indices()[index], frame(), *m_data);
    }

    value_type by_index(size_t const index)
    {
        if (!contains(index))
        {
            throw mol::Error("Atom index " + std::to_string(index) + " not found in the selection");
        }
        return value_type(index, frame(), *m_data);
    }

protected:
    MolData& data()
    {
        return *m_data;
    };

    MolData const& data() const
    {
        return *m_data;
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

        Iterator(MolData& data, indices_iterator begin, Frame frame)
        : m_frame(frame)
        , m_data(&data)
        , m_current(begin)
        {}

        Iterator(Iterator&&) = default;
        Iterator(Iterator const& other) = default;

        value_type operator*() const
        {
            return ItType(*m_current, m_frame, *m_data);
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

        Iterator& operator--()
        {
            m_current--;
            return *this;
        }

        Iterator operator--(int)
        {
            Iterator tmp = *this;
            --(*this);
            return tmp;
        }

        difference_type operator-(const Iterator& other)
        {
            return m_current - other.m_current;
        }

        Iterator& operator+=(difference_type n)
        {
            m_current += n;
            return *this;
        }

        Iterator& operator-=(difference_type n)
        {
            m_current -= n;
            return *this;
        }

        bool operator==(const Iterator& other) const
        {
            return m_current == other.m_current;
        };

        bool operator!=(const Iterator& other) const
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
