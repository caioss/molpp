#include "ResidueBondGuesser.hpp"
#include "tables/ResiduesTable.hpp"
#include <molpp/Bond.hpp>
#include <molpp/AtomSel.hpp>
#include <molpp/ResidueSel.hpp>
#include <algorithm>

using namespace mol::internal;

void ResidueBondGuesser::apply(std::shared_ptr<ResidueSel> residues) const
{
    std::vector<int> bonds_map(RESIDUES_TABLE.max_atoms());

    for (Residue res : *residues)
    {
        if (!RESIDUES_TABLE.contains(res.resname()))
        {
            continue;
        }

        auto const &res_info = RESIDUES_TABLE[res.resname()];
        std::fill(bonds_map.begin(), bonds_map.end(), -1);

        auto atoms = res.atoms();
        for (size_t i = 0; i < (*atoms).size(); i++)
        {
            Atom const atom = (*atoms)[i];
            int const atom_index = res_info.atom_index(atom.name());
            if (atom_index >= 0)
            {
                bonds_map[atom_index] = i;
            }
        }

        for (auto const &bond_info : res_info.bonds)
        {
            int const atom1 = bonds_map[bond_info.atom1];
            int const atom2 = bonds_map[bond_info.atom2];

            if (atom1 < 0 || atom2 < 0)
            {
                // Bond's atoms not present
                continue;
            }

            auto bond = (*atoms)[atom1].bond((*atoms)[atom2]);
            if (!bond)
            {
                // Add a guessed bond
                bond = (*atoms)[atom1].add_bond((*atoms)[atom2]);
                bond->set_guessed(true);
                bond->set_order(bond_info.order);
                bond->set_guessed_order(true);
            }
            else
            {
                // Just fill the missing parameters
                if (bond->order() <= 0)
                {
                    bond->set_order(bond_info.order);
                    bond->set_guessed_order(true);
                };
            }
            bond->set_aromatic(bond_info.aromatic);
        }
    }
}
