#ifndef MOLPP_ATOMSEL_HPP
#define MOLPP_ATOMSEL_HPP

#include <molpp/internal/Sel.hpp>
#include <molpp/Atom.hpp>

#include <vector>

namespace mol
{

class Bond;

class AtomSel : public internal::Sel<Atom, AtomSel>
{
public:
    AtomSel() = delete;
    using internal::Sel<Atom, AtomSel>::Sel;

    auto positions()
    {
        return data()->trajectory().timestep(*frame()).coords()(Eigen::all, indices());
    }

    AtomSel bonded();
    std::vector<std::shared_ptr<mol::Bond>> bonds();

    indices_type const& as_atom_indices() const;

    static indices_type from_atom_indices(internal::IndexRange auto atom_indices, internal::MolData const& /*data*/)
    {
        return atom_indices;
    }

    template<class, class>
    friend class internal::Sel;
};

} // namespace mol

#endif // MOLPP_ATOMSEL_HPP
