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

    auto const as_atom_indices() const
    {
        return data().residues().atom_indices(index());
    }
};

} // namespace mol

#endif // MOLPP_RESIDUE_HPP
