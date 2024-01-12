#ifndef ATOMAGGREGATE_HPP
#define ATOMAGGREGATE_HPP

#include <molpp/Property.hpp>
#include <molpp/internal/requirements.hpp>
#include <molpp/internal/MolData.hpp>

#include <memory>
#include <vector>
#include <concepts>

namespace mol::internal {

template <class Derived>
class AtomAggregate
{
public:
    using coords_type = Eigen::IndexedView<Coord3, Eigen::internal::AllRange<3>, std::vector<index_t>>;
    using const_coords_type = const Eigen::IndexedView<Coord3, Eigen::internal::AllRange<3>, std::vector<index_t>>;

    AtomAggregate() = default;

    AtomAggregate(index_t const index, Frame const frame, internal::MolData* data)
    requires IsAtomAggregate<Derived>
    : m_index{index}
    , m_frame(frame)
    , m_data{data}
    {}

    bool operator==(AtomAggregate<Derived> const &other) const
    {
        return m_data == other.m_data
               && m_index == other.m_index
               && m_frame == other.m_frame;
    }

    //! Index is always read-only
    index_t index() const
    {
        return m_index;
    }

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

    std::vector<std::shared_ptr<Bond>> bonds()
    {
        Derived &derived = static_cast<Derived &>(*this);
        std::vector<index_t> const& indices = derived.atom_indices();
        return m_data->bonds().bonds(indices.begin(), indices.end());
    }

    operator bool() const
    {
        Derived const& derived = static_cast<Derived const&>(*this);
        return m_data && derived.validate_index();
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
        return m_data->properties().template get<Derived, PropertyType>(frame());
    }

    template <IsProperty PropertyType>
    PropertyType const* property() const
    {
        return m_data->properties().template get<Derived, PropertyType>(frame());
    }

    template <IsProperty PropertyType>
    PropertyType::reference get()
    {
        return property<PropertyType>()->value(index());
    }

    template <IsProperty PropertyType>
    PropertyType::const_reference get() const
    {
        return property<PropertyType>()->value(index());
    }

    template <IsProperty PropertyType>
    PropertyType::reference get(PropertyType* attribute)
    {
        return attribute->value(index());
    }

    template <IsProperty PropertyType>
    PropertyType::const_reference get(PropertyType const* attribute) const
    {
        return attribute->value(index());
    }

    template <IsProperty PropertyType>
    void set(PropertyType::value_type const& value)
    {
        get<PropertyType>() = value;
    }

    template <IsProperty PropertyType>
    void set(PropertyType* attribute, PropertyType::value_type const& value)
    {
        get(attribute) = value;
    }

protected:
    MolData* data()
    {
        return m_data;
    }

    MolData const* data() const
    {
        return m_data;
    }

private:
    index_t m_index;
    Frame m_frame;
    internal::MolData* m_data;

    template <class, class>
    friend class Sel;
};

} // namespace mol::internal

#endif // ATOMAGGREGATE_HPP
