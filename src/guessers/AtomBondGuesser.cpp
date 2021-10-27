#include "guessers/AtomBondGuesser.hpp"
#include "tools/SpatialSearch.hpp"
#include <molpp/Bond.hpp>
#include <molpp/AtomSel.hpp>
#include <molpp/MolppCore.hpp>
#include <molpp/ElementsTable.hpp>

using namespace mol::internal;

template <class T>
T pow2(T value)
{
    return value * value;
}

void AtomBondGuesser::apply(AtomSel &atoms) const
{
    auto const coords = atoms.coords();
    float const max_bond_length = 3.0;
    SpatialSearch<AtomSel::coords_type> search(coords, max_bond_length + 0.1);

    for (auto &[atom1, atom2, distance_sq] : search.pairs(max_bond_length))
    {
        int const atomic1 = atoms[atom1].atomic();
        int const atomic2 = atoms[atom2].atomic();

        if (atomic1 == 0 || atomic2 == 0)
        {
            // Unknown atom types
            continue;
        }

        float const radius1 = ELEMENTS_TABLE.covalent_radius(atomic1);
        float const radius2 = ELEMENTS_TABLE.covalent_radius(atomic2);

        // Rule taken from Zhang et al (DOI: 10.1186/1758-2946-4-26)
        float const coff_sq = pow2(radius1 + radius2 + 0.4);
        if (distance_sq > 0.16 && distance_sq < coff_sq)
        {
            auto bond = atoms[atom1].bond(atoms[atom2]);
            if (!bond)
            {
                // Add a guessed bond
                bond = atoms[atom1].add_bond(atoms[atom2]);
                bond->set_guessed(true);
                bond->set_order(1);
                bond->set_guessed_order(true);
            }
        }
    }
}
