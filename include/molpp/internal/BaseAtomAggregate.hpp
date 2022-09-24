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

    BaseAtomAggregate() = delete;
    BaseAtomAggregate(index_t const index, Frame const frame, std::shared_ptr<internal::MolData> data)
    : m_index { index },
      m_frame(frame),
      m_data(data)
    {}

    bool operator==(BaseAtomAggregate const &other) const
    {
        return m_data == other.m_data
               && m_index == other.m_index
               && m_frame == other.m_frame;
    }

    index_t index() const
    {
        return m_index;
    }

    // Read-only
    Frame frame() const
    {
        return m_frame;
    }

protected:
    coords_type coords(std::vector<index_t> &&atom_indices);
    std::vector<std::shared_ptr<Bond>> bonds(std::vector<index_t> const &atom_indices);

    std::shared_ptr<internal::MolData> data()
    {
        return m_data;
    }

    std::shared_ptr<internal::MolData> const data() const
    {
        return m_data;
    }

private:
    index_t m_index;
    Frame m_frame;
    std::shared_ptr<internal::MolData> m_data;
};


} // namespace internal
} // namespace mol

#endif // BASEATOMAGGREGATE_HPP
