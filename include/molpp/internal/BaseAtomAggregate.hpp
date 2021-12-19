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
    using coords_type = Eigen::IndexedView<Coord3, Eigen::internal::AllRange<3>, std::vector<size_t>>;

    BaseAtomAggregate() = delete;
    BaseAtomAggregate(size_t const index, std::optional<size_t> const frame, std::shared_ptr<internal::MolData> data)
    : m_index { index },
      m_frame(frame),
      m_data(data)
    {}

    bool operator==(BaseAtomAggregate const &other) const;

    size_t index() const
    {
        return m_index;
    }

    // Read-only
    std::optional<size_t> frame() const
    {
        return m_frame;
    }

protected:
    coords_type coords(std::vector<size_t> &&atom_indices);
    std::vector<std::shared_ptr<Bond>> bonds(std::vector<size_t> const &atom_indices);

    std::shared_ptr<internal::MolData> data()
    {
        return m_data;
    }

    std::shared_ptr<internal::MolData> const cdata() const
    {
        return m_data;
    }

private:
    size_t m_index;
    std::optional<size_t> m_frame;
    std::shared_ptr<internal::MolData> m_data;
};


} // namespace internal
} // namespace mol

#endif // BASEATOMAGGREGATE_HPP
