#include "guessers/AtomBondGuesser.hpp"
#include "tools/SpatialSearch.hpp"
#include <molpp/Bond.hpp>
#include <molpp/AtomSel.hpp>
#include <molpp/Common.hpp>
#include <molpp/ElementsTable.hpp>

using namespace mol::internal;

template<class T>
T pow2(T value)
{
    return value * value;
}

void AtomBondGuesser::apply(AtomSel& atoms) const
{
    float const max_bond_length = 3.0;
    auto const coords = atoms.positions();
    SpatialSearch<decltype(coords)> search(coords, max_bond_length + 0.1);
    ElementsTable const& elements_table = ELEMENTS_TABLE();

    for (auto& [atom1, atom2, distance_sq] : search.pairs(max_bond_length))
    {
        int const atomic1 = atoms[atom1].atomic_number();
        int const atomic2 = atoms[atom2].atomic_number();

        if (atomic1 == 0 || atomic2 == 0)
        {
            // Unknown atom types
            continue;
        }

        float const radius1 = elements_table.covalent_radius(atomic1);
        float const radius2 = elements_table.covalent_radius(atomic2);

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
