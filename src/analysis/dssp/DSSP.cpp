// TODO revise the includes
#include "analysis/dssp/DSSP.hpp"

#include <molpp/AtomSel.hpp>
#include <molpp/ResidueSel.hpp>
#include <optional>
#include "structure.hpp"
#include <algorithm>
#include "DSSP.hpp"

mol::SSResidue::~SSResidue()
{
}

mol::SSResidue::SSResidue(Residue& residue)
: m_is_proline{residue.resname() == "PRO"}
, N()
, CA()
, C()
, O()
, m_chain{residue.chain()}
, m_residue{std::make_unique<dssp::MResidue>(m_chain, m_is_proline)}
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

void mol::SSResidue::update(dssp::MResidue* previous)
{
    if (!is_amino_acid())
    {
        return;
    }

    m_residue->index = previous ? previous->index + 1 : 0;
    m_residue->previous = previous;
    m_residue->update_positions(N, CA, C, O);

    if (previous)
    {
        previous->next = m_residue.get();
    }

    // Check for chain breaks
    if (previous && !previous->is_valid_distance(m_residue.get()))
    {
        m_residue->index++;
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

void mol::SSResidue::set_structure(SecondaryStructure const structure)
{
    m_residue->structure = structure;
}

bool mol::SSResidue::is_amino_acid() const
{
    return CA && C && O && N;
}

dssp::MResidue& mol::SSResidue::residue() // TODO needed?
{
    return *m_residue;
}

mol::SecondaryStructure mol::SSResidue::secondary_structure() const
{
    // TODO use Unknown as default for dssp::MResidue when not a protein so the check is not needed
    return is_amino_acid() ? m_residue->structure : mol::Unknown;
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
    dssp::MResidue* previous = nullptr; // TODO move this to MProtein
    for (SSResidue& residue : m_residues)
    {
        residue.set_frame(frame);
        if (residue.is_amino_acid())
        {
            // residue.set_structure(mol::Loop);
            residue.update(previous); // TODO move to MProtein
            protein.add_residue(&residue.residue()); // TODO move to MProtein
            previous = &residue.residue(); // TODO parenthesis needed? Look up
        }
        else
        {
            residue.set_structure(mol::Unknown);
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
