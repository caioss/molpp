#ifndef BASESEL_HPP
#define BASESEL_HPP

#include <molpp/MolppCore.hpp>
#include <molpp/internal/SelIndex.hpp>
#include <memory>
#include <optional>

namespace mol {

class Bond;

namespace internal {

class MolData;

class BaseSel
{
public:
    using coords_type = Eigen::IndexedView<Coord3, Eigen::internal::AllRange<3>, std::vector<size_t>>;

    BaseSel() = delete;
    BaseSel(SelIndex&& indices, std::shared_ptr<MolData> data);

    std::optional<size_t> frame() const
    {
        return m_frame;
    }

    void set_frame(std::optional<size_t> frame);

    size_t size() const
    {
        return m_index.size();
    }

    bool contains(size_t const index) const
    {
        return m_index.contains(index);
    }

    std::vector<size_t> const &indices() const
    {
        return m_index.indices();
    }

protected:
    coords_type coords(std::vector<size_t> &&atom_indices);
    std::vector<size_t> bonded(std::vector<size_t> const &atom_indices) const;
    std::vector<std::shared_ptr<mol::Bond>> bonds(std::vector<size_t> const &atom_indices);

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

    std::shared_ptr<MolData> data()
    {
        return m_data;
    };

    std::shared_ptr<MolData> const data() const
    {
        return m_data;
    };

private:
    void init_frame();

    std::optional<size_t> m_frame;
    std::shared_ptr<MolData> m_data;
    SelIndex m_index;
};

} // namespace internal
} // namespace mol

#endif // BASESEL_HPP
