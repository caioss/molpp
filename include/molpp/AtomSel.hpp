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
    using position_type = Eigen::IndexedView<Coord3, Eigen::internal::AllRange<3>, indices_type>;
    using const_position_type = Eigen::IndexedView<Coord3 const, Eigen::internal::AllRange<3>, indices_type>;

    AtomSel() = delete;
    using internal::Sel<Atom, AtomSel>::Sel;

    position_type positions();
    const_position_type positions() const;

    AtomSel bonded();
    std::vector<std::shared_ptr<mol::Bond>> bonds();

    indices_type const& as_atom_indices() const;

    static auto from_atom_indices(internal::IndexRange auto&& atom_indices, internal::MolData const&)
    {
        return atom_indices;
    }

    template<class, class>
    friend class internal::Sel;
};

} // namespace mol

#endif // MOLPP_ATOMSEL_HPP
