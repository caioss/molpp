#ifndef MOLPP_INTERNAL_MOLECULARENTITY_HPP
#define MOLPP_INTERNAL_MOLECULARENTITY_HPP

#include <molpp/MolppCore.hpp>
#include <molpp/internal/MolData.hpp>

#include <memory>
#include <vector>
#include <concepts>

namespace mol::internal
{

class MolData;

class MolecularEntity
{
public:
    MolecularEntity() = default;
    MolecularEntity(index_t const index, Frame const frame, internal::MolData* data);
    bool operator==(MolecularEntity const& other) const;
    //! Index is always read-only
    index_t index() const;
    Frame frame() const;
    void set_frame(Frame const frame);
    bool is_valid() const;

protected:
    MolData* data();
    MolData const* data() const;

private:
    index_t m_index;
    Frame m_frame;
    internal::MolData* m_data;

    template<class>
    friend class Sel;
};

} // namespace mol::internal

#endif // MOLPP_INTERNAL_MOLECULARENTITY_HPP
