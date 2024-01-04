#ifndef MOLPP_PROPERTY_HPP
#define MOLPP_PROPERTY_HPP

#include <molpp/MolppCore.hpp>

#include <vector>
#include <concepts>

namespace mol
{

namespace internal
{
class PropertyTrajectory;
} // namespace internal

class Property;

template<class T>
concept IsProperty = requires(T t)
{
    std::derived_from<T, Property>;
    {t.value(std::declval<index_t>())} -> std::common_reference_with<typename T::value_type&>;
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

protected:
    virtual void resize(size_t const size) = 0;

    friend class mol::internal::PropertyTrajectory;
};

class Position : public Property
{
public:
    using type = Coord3;
    using value_type = Eigen::Block<Coord3, 3, 1, true>;
    using const_value_type = Eigen::Block<const Coord3, 3, 1, true>;

    value_type value(index_t const index)
    {
        return m_positions.col(index);
    };

    const_value_type value(index_t const index) const
    {
        return m_positions.col(index);
    };

    size_t size() const override
    {
        return m_positions.cols();
    }

    Coord3& positions()
    {
        return m_positions;
    }

    Coord3 const& positions() const
    {
        return m_positions;
    }

protected:
    void resize(size_t const size) override
    {
        m_positions.resize(Eigen::NoChange, size);
    }

private:
    Coord3 m_positions;
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

    value_type const& value(index_t const index) const
    {
        return m_values[index];
    };

    size_t size() const override
    {
        return m_values.size();
    }

protected:
    void resize(size_t const size) override
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
