#ifndef MOLPP_PROPERTYTRAJECTORY_HPP
#define MOLPP_PROPERTYTRAJECTORY_HPP

#include <molpp/Property.hpp>

#include <list>
#include <memory>
#include <functional>

namespace mol::internal
{

class PropertyTrajectory
{
public:
    using make_property_fn = std::function<std::unique_ptr<Property>()>;

    PropertyTrajectory(bool const is_time_based, make_property_fn const make_property);
    size_t is_time_based() const;
    void resize(size_t const size);
    size_t num_frames() const;
    Property* add_frame();
    void remove_frame(Frame const frame);
    Property* get(Frame const frame);

private:
    make_property_fn m_make_property;

    Property* emplace_frame();

    bool m_time_based;
    std::list<std::unique_ptr<Property>> m_trajectory;
};

} // namespace mol::internal

#endif // MOLPP_PROPERTYTRAJECTORY_HPP
