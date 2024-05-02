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

protected:
    static size_t data_size(internal::MolData const& data);
    std::vector<index_t> const& atom_indices() const;

    static internal::SelIndex from_atom_indices(internal::IndexRange auto const& atom_indices, internal::MolData const& data)
    {
        return {atom_indices, data_size(data)};
    }

    template<class, class>
    friend class internal::Sel;
};

} // namespace mol

#endif // MOLPP_ATOMSEL_HPP
