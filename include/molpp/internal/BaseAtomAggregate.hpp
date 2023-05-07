#ifndef BASEATOMAGGREGATE_HPP
#define BASEATOMAGGREGATE_HPP

#include <molpp/MolppCore.hpp>
#include <memory>
#include <vector>
#include <optional>

namespace mol {

class Bond;

namespace internal {

class MolData;

class BaseAtomAggregate
{
public:
    using coords_type = Eigen::IndexedView<Coord3, Eigen::internal::AllRange<3>, std::vector<index_t>>;
    using const_coords_type = const Eigen::IndexedView<Coord3, Eigen::internal::AllRange<3>, std::vector<index_t>>;

    BaseAtomAggregate() = default;

    BaseAtomAggregate(index_t const index, Frame const frame, internal::MolData* data)
    : m_index{index},
      m_frame(frame),
      m_data{data}
    {}

    bool operator==(BaseAtomAggregate const &other) const
    {
        return m_data == other.m_data
               && m_index == other.m_index
               && m_frame == other.m_frame;
    }

    // Read-only
    index_t index() const
    {
        return m_index;
    }

    Frame frame() const;
    void set_frame(Frame const frame);

protected:
    coords_type coords(std::vector<index_t> &&atom_indices);
    const_coords_type coords(std::vector<index_t> &&atom_indices) const;
    std::vector<std::shared_ptr<Bond>> bonds(std::vector<index_t> const &atom_indices);

    internal::MolData* data()
    {
        return m_data;
    }

    internal::MolData const* data() const
    {
        return m_data;
    }

private:
    index_t m_index;
    Frame m_frame;
    internal::MolData* m_data;
};


} // namespace internal
} // namespace mol

#endif // BASEATOMAGGREGATE_HPP
