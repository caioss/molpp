#ifndef ATOM_HPP
#define ATOM_HPP

#include <molpp/MolppCore.hpp>
#include <molpp/internal/AtomAggregate.hpp>
#include <memory>
#include <vector>
#include <optional>

namespace mol {

class Bond;
class Residue;

class Atom : public internal::AtomAggregate<Atom>
{
    friend class internal::AtomAggregate<Atom>;

public:
    using internal::AtomAggregate<Atom>::AtomAggregate;

    int resid() const;
    Residue residue();
    index_t residue_id() const;

    std::shared_ptr<Bond> add_bond(index_t const bonded_to);
    std::shared_ptr<Bond> add_bond(Atom const &bonded_to);
    std::shared_ptr<Bond> bond(index_t const other);
    std::shared_ptr<Bond> bond(Atom const &other);

    std::vector<index_t> atom_indices() const;

protected:
    bool validate_index() const;
};

} // namespace mol

#endif // ATOM_HPP
