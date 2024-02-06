#ifndef MOLPP_PROPERTYCONTAINER_HPP
#define MOLPP_PROPERTYCONTAINER_HPP

#include <molpp/internal/PropertyTrajectory.hpp>
#include <molpp/Property.hpp>
#include <molpp/internal/requirements.hpp>

#include <map>
#include <memory>
#include <utility>
#include <typeinfo>
#include <typeindex>
#include <unordered_map>

namespace mol::internal
{

class PropertyContainer
{
public:
    using property_key_type = std::pair<std::type_index, std::type_index>;
    using size_key_type = std::type_index;

    PropertyContainer();

    size_t num_frames() const;
    Frame add_frame();
    void remove_frame(Frame const frame);

    template<IsAtomAggregate AggregateType>
    void resize(size_t const size);

    template<IsAtomAggregate AggregateType, IsProperty PropertyType>
    PropertyType const* get(Frame const frame) const;

    template<IsAtomAggregate AggregateType, IsProperty PropertyType>
    PropertyType* get(Frame const frame);

    template<IsAtomAggregate AggregateType, IsProperty PropertyType>
    PropertyType* add(bool const is_time_based, size_t const size);

private:
    void resize(size_key_type const key, size_t const size);
    Property* emplace(property_key_type const key, bool const is_time_based, size_t const size, PropertyTrajectory::make_property_fn const make_property);
    Property const* get(property_key_type const key, Frame const frame) const;
    Property* get(property_key_type const key, Frame const frame);

    template<IsAtomAggregate AggregateType, IsProperty PropertyType>
    static std::unique_ptr<Property> make_property();

    size_t m_num_frames;
    std::map<property_key_type, PropertyTrajectory> m_properties;
};

template<IsAtomAggregate AggregateType>
void PropertyContainer::resize(size_t const size)
{
    resize(typeid(AggregateType), size);
}

template<IsAtomAggregate AggregateType, IsProperty PropertyType>
PropertyType const* PropertyContainer::get(Frame const frame) const
{
    property_key_type key(typeid(AggregateType), typeid(PropertyType));
    Property const* property = get(key, frame);

    if (!property)
    {
        return nullptr;
    }
    return static_cast<PropertyType const*>(property);
}

template<IsAtomAggregate AggregateType, IsProperty PropertyType>
PropertyType* PropertyContainer::get(Frame const frame)
{
    // Delegate to the const version
    return const_cast<PropertyType*>(std::as_const(*this).get<AggregateType, PropertyType>(frame));
}

template<IsAtomAggregate AggregateType, IsProperty PropertyType>
PropertyType* PropertyContainer::add(bool const is_time_based, size_t const size)
{
    property_key_type key(typeid(AggregateType), typeid(PropertyType));
    Property* property = emplace(key, is_time_based, size, make_property<AggregateType, PropertyType>);
    return static_cast<PropertyType*>(property);
}

template<IsAtomAggregate AggregateType, IsProperty PropertyType>
std::unique_ptr<Property> PropertyContainer::make_property()
{
    return std::make_unique<PropertyType>();
}

} // namespace mol::internal

#endif // MOLPP_PROPERTYCONTAINER_HPP
