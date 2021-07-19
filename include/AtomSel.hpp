#ifndef ATOMSEL_HPP
#define ATOMSEL_HPP

#include "Atom.hpp"
#include "MolppCore.hpp"
#include "internal/BaseSel.hpp"
#include <vector>
#include <memory>

namespace mol {

class AtomSel : public internal::BaseSel<Atom>
{

public:
    using coords_type = Eigen::IndexedView<Coord3, Eigen::internal::AllRange<3>, std::vector<size_t>>;

    AtomSel() = delete;
    using internal::BaseSel<Atom>::BaseSel;

    coords_type coords();
    std::shared_ptr<AtomSel> bonded();
    std::vector<std::shared_ptr<Bond>> bonds();
};

} // namespace mol

#endif // ATOMSEL_HPP
