#ifndef RESIDUESEL_HPP
#define RESIDUESEL_HPP

#include <molpp/internal/Sel.hpp>
#include <molpp/Residue.hpp>
#include <memory>

namespace mol {

class ResidueSel : public internal::Sel<Residue, ResidueSel>
{
public:
    ResidueSel() = delete;
    using internal::Sel<Residue, ResidueSel>::Sel;
    std::vector<index_t> atom_indices() const;

protected:
    static std::vector<index_t> from_atom_indices(std::vector<index_t> &&atom_indices, std::shared_ptr<mol::internal::MolData> const data);
    static size_t max_size(std::shared_ptr<mol::internal::MolData> const data);

    template <class, class>
    friend class internal::Sel;
};

} // namespace mol

#endif // RESIDUESEL_HPP
