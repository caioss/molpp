#ifndef MOLPP_PROPERTYCONTAINER_HPP
#define MOLPP_PROPERTYCONTAINER_HPP

#include <molpp/internal/PropertyTrajectory.hpp>
#include <molpp/Property.hpp>
#include <molpp/internal/requirements.hpp>

#include <map>
#include <memory>
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
    size_t size() const;

    template<IsAtomAggregate AggregateType>
    void set_size(size_t const size);

    template<IsAtomAggregate AggregateType, IsProperty PropertyType>
    PropertyType* get(Frame const frame);

    template<IsAtomAggregate AggregateType, IsProperty PropertyType>
    PropertyType* add(bool const is_time_based);

private:
    size_t size(size_key_type const key) const;
    void resize(size_key_type const key, size_t const size);
    Property* emplace(property_key_type const key, bool const is_time_based, PropertyTrajectory::make_property_fn const make_property);
    Property* get(property_key_type const key, Frame const frame);

    template<IsAtomAggregate AggregateType, IsProperty PropertyType>
    static std::unique_ptr<Property> make_property();

    size_t m_num_frames;
    std::unordered_map<size_key_type, size_t> m_sizes;
    std::map<property_key_type, PropertyTrajectory> m_properties;
};

template<IsAtomAggregate AggregateType>
size_t PropertyContainer::size() const
{
    return size(typeid(AggregateType));
}

template<IsAtomAggregate AggregateType>
void PropertyContainer::set_size(size_t const size)
{
    resize(typeid(AggregateType), size);
}

template<IsAtomAggregate AggregateType, IsProperty PropertyType>
PropertyType* PropertyContainer::get(Frame const frame)
{
    property_key_type key(typeid(AggregateType), typeid(PropertyType));
    Property* property = get(key, frame);

    if (!property)
    {
        return nullptr;
    }
    return static_cast<PropertyType*>(property);
}

template<IsAtomAggregate AggregateType, IsProperty PropertyType>
PropertyType* PropertyContainer::add(bool const is_time_based)
{
    property_key_type key(typeid(AggregateType), typeid(PropertyType));
    Property* property = emplace(key, is_time_based, make_property<AggregateType, PropertyType>);
    return static_cast<PropertyType*>(property);
}

template<IsAtomAggregate AggregateType, IsProperty PropertyType>
std::unique_ptr<Property> PropertyContainer::make_property()
{
    return std::make_unique<PropertyType>();
}

} // namespace mol::internal

#endif // MOLPP_PROPERTYCONTAINER_HPP
