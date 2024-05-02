#ifndef MOLPP_INTERNAL_ATOMAGGREGATE_HPP
#define MOLPP_INTERNAL_ATOMAGGREGATE_HPP

#include <molpp/internal/requirements.hpp>
#include <memory>
#include <vector>
#include <concepts>
#include <molpp/internal/MolData.hpp>

namespace mol::internal
{

class AtomAggregate
{
public:
    AtomAggregate() = default;

    AtomAggregate(index_t const index, Frame const frame, internal::MolData* data)
    : m_index{index}
    , m_frame(frame)
    , m_data{data}
    {}

    bool operator==(AtomAggregate const& other) const
    {
        return m_data == other.m_data && m_index == other.m_index && m_frame == other.m_frame;
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
        if (frame && frame >= m_data->trajectory().num_frames())
        {
            throw mol::MolError("Out of bounds frame: " + std::to_string(*frame));
        }
        m_frame = frame;
    }

    bool is_valid() const
    {
        return m_data;
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

    template<class, class>
    friend class Sel;
};

} // namespace mol::internal

#endif // MOLPP_INTERNAL_ATOMAGGREGATE_HPP
