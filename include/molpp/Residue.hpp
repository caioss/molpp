#ifndef MOLPP_RESIDUE_HPP
#define MOLPP_RESIDUE_HPP

#include <molpp/MolppCore.hpp>
#include <molpp/internal/MolecularEntity.hpp>

#include <optional>

namespace mol
{

class Atom;
class Chain;

class Residue : public internal::MolecularEntity
{
public:
    using internal::MolecularEntity::MolecularEntity;

    static MolecularEntityCategory category();

    std::optional<Chain> chain();
    std::optional<index_t> chain_index() const;

    int id() const;
    void set_id(int const resid);

    std::string const& name() const;
    void set_name(std::string const& resname);

    std::string const& segid() const;
    void set_segid(std::string const& segid);

    void add_atom(index_t index);
    void add_atom(Atom const& atom);
    size_t size() const;
};

} // namespace mol

#endif // MOLPP_RESIDUE_HPP
