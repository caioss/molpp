#ifndef ATOMSEL_HPP
#define ATOMSEL_HPP

#include <molpp/internal/Sel.hpp>
#include <molpp/Atom.hpp>
#include <vector>
#include <memory>

namespace mol {

class Bond;

class AtomSel : public internal::Sel<Atom, AtomSel>
{
public:
    AtomSel() = delete;
    using internal::Sel<Atom, AtomSel>::Sel;
    std::vector<index_t> atom_indices() const;

protected:
    static std::vector<index_t> from_atom_indices(std::vector<index_t> &&atom_indices, std::shared_ptr<mol::internal::MolData> const data);
    static size_t max_size(std::shared_ptr<mol::internal::MolData> const data);

    template <class, class>
    friend class internal::Sel;
};

} // namespace mol

#endif // ATOMSEL_HPP
