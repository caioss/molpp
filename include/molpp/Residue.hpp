#ifndef MOLPP_RESIDUE_HPP
#define MOLPP_RESIDUE_HPP

#include <molpp/MolppCore.hpp>
#include <molpp/internal/MolecularEntity.hpp>

namespace mol
{

class Atom;

class Residue : public internal::MolecularEntity
{
public:
    using internal::MolecularEntity::MolecularEntity;

    static MolecularEntityCategory category();

    int id() const;
    void set_id(int const resid);

    std::string const& name() const;
    void set_name(std::string const& resname);

    std::string const& segid() const;
    void set_segid(std::string const& segid);

    std::string const& chain() const;
    void set_chain(std::string const& chain);

    void add_atom(index_t index);
    void add_atom(Atom const& atom);
    size_t size() const;
};

} // namespace mol

#endif // MOLPP_RESIDUE_HPP
