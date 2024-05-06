#ifndef MOLPP_RESIDUE_HPP
#define MOLPP_RESIDUE_HPP

#include <molpp/MolppCore.hpp>
#include <molpp/internal/AtomAggregate.hpp>

namespace mol
{

class Atom;

class Residue : public internal::AtomAggregate
{
public:
    using internal::AtomAggregate::AtomAggregate;

    int residue_id() const;
    void set_residue_id(int const resid);

    std::string const& residue_name() const;
    void set_residue_name(std::string const& resname);

    std::string const& segid() const;
    void set_segid(std::string const& segid);

    std::string const& chain() const;
    void set_chain(std::string const& chain);

    void add_atom(index_t index);
    void add_atom(Atom const& atom);
    size_t size() const;

    std::vector<index_t> as_atom_indices() const;
};

} // namespace mol

#endif // MOLPP_RESIDUE_HPP
