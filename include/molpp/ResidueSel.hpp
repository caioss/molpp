#ifndef RESIDUESEL_HPP
#define RESIDUESEL_HPP

#include <molpp/internal/Sel.hpp>
#include <molpp/Residue.hpp>
#include <unordered_set>

namespace mol {

class ResidueSel : public internal::Sel<Residue, ResidueSel>
{
public:
    ResidueSel() = delete;
    using internal::Sel<Residue, ResidueSel>::Sel;

protected:
    static size_t data_size(internal::MolData const& data);
    std::vector<index_t> atom_indices() const;

    static internal::SelIndex from_atom_indices(internal::IndexRange auto const& atom_indices, internal::MolData const& data)
    {
        std::unordered_set<index_t> residues;
        for (auto const index : atom_indices)
        {
            residues.insert(atom_index(index, data));
        }

        return internal::SelIndex(residues, data_size(data));
    }

private:
    static size_t atom_index(size_t const atom_index, internal::MolData const& data);

    template <class, class>
    friend class internal::Sel;
};

} // namespace mol

#endif // RESIDUESEL_HPP
