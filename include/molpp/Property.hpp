#ifndef MOLPP_PROPERTY_HPP
#define MOLPP_PROPERTY_HPP

#include <molpp/MolppCore.hpp>

#include <vector>

namespace mol
{

class Property;

template<class T>
concept IsProperty = requires(T t)
{
    std::derived_from<T, Property>;
    {t.value(std::declval<index_t>())} -> std::same_as<typename T::value_type&>;
};

template<class T>
concept IsContainerProperty = requires(T t)
{
    {t[std::declval<index_t>()]} -> std::same_as<typename T::value_type&>;
    {t.size()} -> std::same_as<size_t>;
    {t.resize(std::declval<index_t>())} -> std::same_as<void>;
};

class Property
{
public:
    virtual ~Property() = default;
    virtual size_t size() const = 0;
    virtual void resize(size_t const size) = 0;
};

template<class Type, class Container = std::vector<Type>>
class ContainerProperty : public Property
{
public:
    using value_type = typename Container::value_type;

    value_type& value(index_t const index)
    {
        return m_values[index];
    };

    size_t size() const
    {
        return m_values.size();
    }

    void resize(size_t const size)
    {
        m_values.resize(size);
    }

private:
    Container m_values;
};

class Occupancy : public ContainerProperty<float>{};
class TemperatureFactor : public ContainerProperty<float>{};
class Mass : public ContainerProperty<float>{};
class Charge : public ContainerProperty<float>{};
class Radius : public ContainerProperty<float>{};
class AtomicNumber : public ContainerProperty<int>{};
class ResID : public ContainerProperty<int>{};
class Name : public ContainerProperty<std::string>{};
class Type : public ContainerProperty<std::string>{};
class AlternateLocation : public ContainerProperty<std::string>{};
class InsertionCode : public ContainerProperty<std::string>{};
class ResName : public ContainerProperty<std::string>{};

}; // namespace mol

#endif // MOLPP_PROPERTY_HPP
