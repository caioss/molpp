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

    std::vector<size_t> atom_indices() const;
    static ResidueSel from_atom_indices(std::vector<size_t> &&atom_indices, std::shared_ptr<mol::internal::MolData> data);

protected:
    static size_t max_size(std::shared_ptr<mol::internal::MolData> data);

    friend class internal::Sel<Residue, ResidueSel>;
};

} // namespace mol

#endif // RESIDUESEL_HPP
