#include "ResidueBondGuesser.hpp"
#include "tables/ResiduesTable.hpp"
#include <molpp/Bond.hpp>
#include <molpp/AtomSel.hpp>
#include <molpp/ResidueSel.hpp>

#include <optional>
#include <algorithm>

namespace mol::internal
{

std::optional<int> find_atom(ResiduesTable::Residue const& residue_info, std::string const& name)
{
    auto const atom_it = residue_info.atoms.find(name);
    if (atom_it == residue_info.atoms.end())
    {
        return std::nullopt;
    }
    return atom_it->second;
}

void ResidueBondGuesser::apply(ResidueSel& residue_sel) const
{
    ResiduesTable const& residues_table = RESIDUES_TABLE();
    std::vector<int> bonds_map;

    for (Residue residue : residue_sel)
    {
        if (!residues_table.contains(residue.name()))
        {
            continue;
        }

        ResiduesTable::Residue const& residue_info = residues_table[residue.name()];

        // Clear bonds map
        std::fill(bonds_map.begin(), bonds_map.end(), -1);
        bonds_map.resize(residue_info.atoms.size(), -1);

        AtomSel atoms(residue);
        for (index_t i = 0; i < atoms.size(); i++)
        {
            Atom const atom = atoms[i];
            std::optional<int> const atom_index = find_atom(residue_info, atom.name());
            if (atom_index)
            {
                bonds_map[*atom_index] = i;
            }
        }

        for (auto const& bond_info : residue_info.bonds)
        {
            int const atom1 = bonds_map[bond_info.atom1];
            int const atom2 = bonds_map[bond_info.atom2];

            if (atom1 < 0 || atom2 < 0)
            {
                // Bonded atoms not present
                continue;
            }

            std::shared_ptr<mol::Bond> bond = atoms[atom1].bond(atoms[atom2]);
            if (!bond)
            {
                // Add a guessed bond
                bond = atoms[atom1].add_bond(atoms[atom2]);
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

} // namespace mol::internal
