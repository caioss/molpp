#ifndef MOLPP_CHAIN_HPP
#define MOLPP_CHAIN_HPP

#include <molpp/MolppCore.hpp>
#include <molpp/internal/MolecularEntity.hpp>

namespace mol
{

class Residue;

class Chain : public internal::MolecularEntity
{
public:
    using internal::MolecularEntity::MolecularEntity;

    static MolecularEntityCategory category();
    size_t size() const;

    std::string const& name() const;
    void set_name(std::string const& resname);

    void add_residue(index_t index);
    void add_residue(Residue const& residue);
};

} // namespace mol

#endif // MOLPP_CHAIN_HPP
