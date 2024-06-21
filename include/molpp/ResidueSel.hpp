#ifndef MOLPP_RESIDUESEL_HPP
#define MOLPP_RESIDUESEL_HPP

#include <molpp/internal/Sel.hpp>
#include <molpp/Residue.hpp>

namespace mol
{

class ResidueSel;

namespace internal
{

template<>
struct SelTraits<ResidueSel>
{
    using entity_type = Residue;
};

} // namespace internal

class ResidueSel : public internal::Sel<ResidueSel>
{
public:
    ResidueSel() = delete;
    using internal::Sel<ResidueSel>::Sel;

    indices_type as_atom_indices() const;

    static indices_type from_atom_indices(internal::IndexRange auto const& atom_indices, internal::MolData const& data)
    {
        ResidueSel::indices_type residues;
        for (auto const index : atom_indices)
        {
            internal::Topology const& topology = data.topology();
            std::optional<index_t> const residue_index = topology.find_link({MolecularEntityCategory::Atom, index}, MolecularEntityCategory::Residue);
            if (residue_index)
            {
                residues.push_back(*residue_index);
            }
        }

        return residues;
    }

    template<class>
    friend class internal::Sel;
};

} // namespace mol

#endif // MOLPP_RESIDUESEL_HPP
