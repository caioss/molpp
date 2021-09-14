#ifndef BASESEL_HPP
#define BASESEL_HPP

#include <molpp/MolppCore.hpp>
#include <molpp/MolError.hpp>
#include <vector>
#include <memory>
#include <optional>

namespace mol {

class Bond;
class AtomSel;
class ResidueSel;

namespace internal {

class MolData;

class BaseSel
{
public:
    using coords_type = Eigen::IndexedView<Coord3, Eigen::internal::AllRange<3>, std::vector<size_t>>;

    BaseSel() = delete;
    BaseSel(size_t const max_size, std::shared_ptr<MolData> data);
    BaseSel(size_t const max_size, std::vector<size_t> const &indices, std::shared_ptr<MolData> data);

    std::optional<size_t> frame() const
    {
        return m_frame;
    }

    void set_frame(std::optional<size_t> frame);

    size_t size() const
    {
        return m_indices.size();
    }

    bool contains(size_t const index) const
    {
        if (index >= m_selected.size())
        {
            return false;
        }
        else
        {
            return m_selected[index];
        }
    }

    std::vector<size_t> const &indices() const
    {
        return m_indices;
    }

    std::vector<bool> const &selected() const
    {
        return m_selected;
    }

protected:
    coords_type coords(std::vector<size_t> const &atom_indices);
    std::vector<size_t> bonded(std::vector<size_t> const &atom_indices) const;
    std::vector<std::shared_ptr<mol::Bond>> bonds(std::vector<size_t> const &atom_indices);
    std::shared_ptr<AtomSel> atoms(std::vector<size_t> &&atom_indices);
    std::shared_ptr<ResidueSel> residues(std::vector<size_t> &&atom_indices);

    std::vector<size_t>::iterator indices_begin()
    {
        return m_indices.begin();
    }

    std::vector<size_t>::iterator indices_end()
    {
        return m_indices.end();
    }

    std::shared_ptr<MolData> const cdata() const
    {
        return m_data;
    };

    std::shared_ptr<MolData> data()
    {
        return m_data;
    };

private:
    void init_frame();

    std::optional<size_t> m_frame;
    std::shared_ptr<MolData> m_data;
    std::vector<bool> m_selected;
    std::vector<size_t> m_indices;
};

} // namespace internal
} // namespace mol

#endif // BASESEL_HPP
