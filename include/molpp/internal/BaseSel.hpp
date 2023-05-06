#ifndef BASESEL_HPP
#define BASESEL_HPP

#include <molpp/MolppCore.hpp>
#include <molpp/internal/SelIndex.hpp>
#include <molpp/Timestep.hpp>
#include <memory>
#include <optional>

namespace mol {

class Bond;

namespace internal {

class MolData;

class BaseSel
{
public:
    using coords_type = Eigen::IndexedView<Coord3, Eigen::internal::AllRange<3>, std::vector<index_t>>;

    BaseSel() = delete;
    BaseSel(BaseSel &&) = default;
    BaseSel(BaseSel const&) = default;
    BaseSel& operator=(BaseSel &&) = default;
    BaseSel& operator=(BaseSel const&) = default;
    BaseSel(SelIndex&& indices, MolData* data);

    Frame frame() const
    {
        return m_frame;
    }

    void set_frame(Frame frame);

    size_t size() const
    {
        return m_index.size();
    }

    bool contains(index_t const index) const
    {
        return m_index.contains(index);
    }

    std::vector<index_t> const &indices() const
    {
        return m_index.indices();
    }

    Timestep& timestep();

protected:
    std::vector<index_t> bonded(std::vector<index_t> const &atom_indices) const;
    std::vector<std::shared_ptr<mol::Bond>> bonds(std::vector<index_t> const &atom_indices);

    SelIndex::iterator indices_begin()
    {
        return m_index.indices_begin();
    }

    SelIndex::iterator indices_begin() const
    {
        return m_index.indices_begin();
    }

    SelIndex::iterator indices_end()
    {
        return m_index.indices_end();
    }

    SelIndex::iterator indices_end() const
    {
        return m_index.indices_end();
    }

    MolData* data()
    {
        return m_data;
    };

    MolData const* data() const
    {
        return m_data;
    };

private:
    void init_frame();

    Frame m_frame;
    MolData* m_data;
    SelIndex m_index;
};

} // namespace internal
} // namespace mol

#endif // BASESEL_HPP
