#ifndef RESIDUESTABLE_HPP
#define RESIDUESTABLE_HPP

#include <molpp/MolppCore.hpp>
#include <string>
#include <vector>
#include <unordered_map>
#include <initializer_list>

namespace mol::internal {

class ResiduesTable
{
public:
    struct Bond
    {
        bool aromatic;
        int order;
        index_t atom1;
        index_t atom2;
    };

    struct Residue
    {
        std::unordered_map<std::string, index_t> atoms;
        std::vector<Bond> bonds;

        int atom_index(std::string const &name) const
        {
            auto const atom_it = atoms.find(name);
            if (atom_it == atoms.end())
            {
                return -1;
            }
            return atom_it->second;
        }
    };

    ResiduesTable(std::initializer_list<std::pair<std::string, Residue>> data);

    size_t max_atoms() const
    {
        return m_max_atoms;
    }

    bool contains(std::string const &resname) const
    {
        return m_residues.contains(resname);
    }

    Residue const &operator[](std::string const &resname) const
    {
        return m_residues.at(resname);
    }

private:
    size_t m_max_atoms;
    std::unordered_map<std::string, Residue> m_residues;
};

// Residues data derived from RCSB PDB Ligand Expo:
// Ligand Depot: a data warehouse for ligands bound to macromolecules.
// Bioinformatics. 2004 Sep 1;20(13):2153-5.
// Feng Z, Chen L, Maddula H, Akcan O, Oughtred R, Berman HM, Westbrook J.
// PubMed: 15059838
ResiduesTable const& RESIDUES_TABLE();

} // namespace mol::internal

#endif // RESIDUESTABLE_HPP
