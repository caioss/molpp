#ifndef MOLPP_ATOMSEL_HPP
#define MOLPP_ATOMSEL_HPP

#include <molpp/internal/Sel.hpp>
#include <molpp/Atom.hpp>

#include <vector>

namespace mol
{

class Bond;
class AtomSel;

namespace internal
{

template<>
struct SelTraits<AtomSel>
{
    using entity_type = Atom;
};

} // namespace internal

class AtomSel : public internal::Sel<AtomSel>
{
public:
    using position_type = Eigen::IndexedView<Coord3, Eigen::internal::AllRange<3>, indices_type>;
    using const_position_type = Eigen::IndexedView<Coord3 const, Eigen::internal::AllRange<3>, indices_type>;

    AtomSel() = delete;
    using internal::Sel<AtomSel>::Sel;

    position_type positions();
    const_position_type positions() const;

    AtomSel bonded();
    std::vector<std::shared_ptr<mol::Bond>> bonds();

    template<class>
    friend class internal::Sel;
};

} // namespace mol

#endif // MOLPP_ATOMSEL_HPP
