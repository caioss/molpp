#ifndef RESIDUE_HPP
#define RESIDUE_HPP

#include <molpp/MolppCore.hpp>
#include <molpp/internal/AtomAggregate.hpp>

namespace mol {

class Atom;

class Residue : public internal::AtomAggregate<Residue>
{
public:
    Residue() = delete;
    using internal::AtomAggregate<Residue>::AtomAggregate;

    int resid() const;
    void set_resid(int const &resid);

    std::string resname() const;
    void set_resname(std::string const &resname);

    std::string segid() const;
    void set_segid(std::string const &segid);

    std::string chain() const;
    void set_chain(std::string const &chain);

    void add_atom(index_t index);
    void add_atom(Atom const &atom);
    size_t size() const;

    std::vector<index_t> atom_indices() const;
};

} // namespace mol

#endif // RESIDUE_HPP
