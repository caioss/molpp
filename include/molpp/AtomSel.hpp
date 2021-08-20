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

    std::vector<size_t> atom_indices() const;
    static AtomSel from_atom_indices(std::vector<size_t> &&atom_indices, std::shared_ptr<mol::internal::MolData> data);

protected:
    static size_t max_size(std::shared_ptr<mol::internal::MolData> data);

    friend class internal::Sel<Atom, AtomSel>;
};

} // namespace mol

#endif // ATOMSEL_HPP
