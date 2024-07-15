#ifndef MOLPP_TABLES_RESIDUESTABLE_HPP
#define MOLPP_TABLES_RESIDUESTABLE_HPP

#include <molpp/Common.hpp>

#include <string>
#include <vector>
#include <unordered_map>
#include <initializer_list>

namespace mol::internal
{

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
    };

    ResiduesTable(std::initializer_list<std::pair<std::string, Residue>> data);

    bool contains(std::string const& resname) const
    {
        return m_residues.contains(resname);
    }

    Residue const& operator[](std::string const& resname) const
    {
        return m_residues.at(resname);
    }

private:
    std::unordered_map<std::string, Residue> m_residues;
};

// Residues data derived from RCSB PDB Ligand Expo:
// Ligand Depot: a data warehouse for ligands bound to macromolecules.
// Bioinformatics. 2004 Sep 1;20(13):2153-5.
// Feng Z, Chen L, Maddula H, Akcan O, Oughtred R, Berman HM, Westbrook J.
// PubMed: 15059838
ResiduesTable const& RESIDUES_TABLE();

} // namespace mol::internal

#endif // MOLPP_TABLES_RESIDUESTABLE_HPP
