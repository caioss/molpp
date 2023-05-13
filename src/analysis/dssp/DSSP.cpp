// TODO revise the includes
#include "analysis/dssp/DSSP.hpp"

#include <molpp/AtomSel.hpp>
#include <molpp/ResidueSel.hpp>
#include <optional>
#include "structure.hpp"
#include <algorithm>
#include "DSSP.hpp"

mol::SSResidue::SSResidue(Residue& residue)
: m_is_proline{residue.resname() == "PRO"}
, N()
, CA()
, C()
, O()
, m_chain{residue.chain()}
, m_residue{nullptr}
{
    for (Atom atom : AtomSel(residue))
    {
        std::string const name = atom.name();
        if (name == "N")
        {
            N = atom;
        }
        else if (name == "CA")
        {
            CA = atom;
        }
        else if (name == "C")
        {
            C = atom;
        }
        else if (name == "O")
        {
            O = atom;
        }
    }
}

void mol::SSResidue::set_frame(Frame const frame)
{
    if (N)
    {
        N.set_frame(frame);
    }
    if (CA)
    {
        CA.set_frame(frame);
    }
    if (C)
    {
        C.set_frame(frame);
    }
    if (O)
    {
        O.set_frame(frame);
    }
}

bool mol::SSResidue::is_amino_acid() const
{
    return CA && C && O && N;
}

mol::SecondaryStructure mol::SSResidue::secondary_structure() const
{
    return m_residue ? m_residue->structure : mol::Unknown;
}

bool mol::SSResidue::is_proline() const
{
    return m_is_proline;
}

std::string const& mol::SSResidue::chain() const
{
    return m_chain;
}

mol::DSSP::DSSP(MolSystem& molecule)
{
    ResidueSel residues(molecule.atoms(0));
    for (mol::Residue residue : residues)
    {
        m_residues.emplace_back(residue);
    }
}

std::vector<mol::SecondaryStructure> mol::DSSP::run(mol::Frame frame)
{
    dssp::MProtein protein;
    dssp::MResidue* previous = nullptr;
    for (SSResidue& residue : m_residues)
    {
        residue.set_frame(frame);
        if (residue.is_amino_acid())
        {
            dssp::MResidue* mresidue = protein.emplace_residue(residue.chain(), residue.is_proline(), previous, residue.N.coords(), residue.CA.coords(), residue.C.coords(), residue.O.coords());
            residue.m_residue = mresidue;
            previous = mresidue;
        }
    }

    protein.CalculateSecondaryStructure();

    std::vector<mol::SecondaryStructure> structures;
    structures.reserve(protein.residues().size());

    for (SSResidue const& residue : m_residues)
    {
            structures.push_back(residue.secondary_structure());
    }

    return structures;
}
