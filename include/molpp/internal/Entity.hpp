#ifndef MOLPP_INTERNAL_ENTITY_HPP
#define MOLPP_INTERNAL_ENTITY_HPP

#include <molpp/Common.hpp>
#include <molpp/internal/MolData.hpp>

#include <memory>
#include <vector>
#include <concepts>

namespace mol::internal
{

class MolData;

class Entity
{
public:
    Entity() = delete;
    Entity(index_t const index, Frame const frame, internal::MolData& data);
    bool operator==(Entity const& other) const;
    //! Index is always read-only
    index_t index() const;
    std::vector<index_t> indices() const;
    Frame frame() const;
    void set_frame(Frame const frame);

protected:
    MolData& data();
    MolData const& data() const;

private:
    index_t m_index;
    Frame m_frame;
    MolData* m_data;

    template<class>
    friend class Sel;
};

} // namespace mol::internal

#endif // MOLPP_INTERNAL_ENTITY_HPP
