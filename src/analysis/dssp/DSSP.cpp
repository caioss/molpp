#include "analysis/dssp/DSSP.hpp"
#include "analysis/dssp/structure.hpp"
#include <molpp/AtomSel.hpp>
#include <molpp/ResidueSel.hpp>
#include <molpp/Property.hpp>

mol::SSResidue::SSResidue()
: m_is_proline{false}
, m_is_chain_break{false}
, m_structure{mol::Unknown}
, m_N()
, m_CA()
, m_C()
, m_O()
{
}

mol::SSResidue::SSResidue(Residue& residue)
: m_is_chain_break{false}
, m_structure{mol::Unknown}
, m_N()
, m_CA()
, m_C()
, m_O()
, m_chain{residue.chain()}
{
    AtomSel atom_sel(residue);
    Name* name_property = atom_sel.property<Name>();

    for (Atom atom : atom_sel)
    {
        std::string const name = atom.get(name_property);
        if (name == "N")
        {
            m_N = atom;
        }
        else if (name == "CA")
        {
            m_CA = atom;
        }
        else if (name == "C")
        {
            m_C = atom;
        }
        else if (name == "O")
        {
            m_O = atom;
        }
    }

    if (is_amino_acid())
    {
        m_is_proline = residue.resname() == "PRO";
    }
    else
    {
        m_N = Atom();
        m_CA = Atom();
        m_C = Atom();
        m_O = Atom();
    }
}

void mol::SSResidue::set_frame(Frame const frame)
{
    if (m_N)
    {
        m_N.set_frame(frame);
    }
    if (m_CA)
    {
        m_CA.set_frame(frame);
    }
    if (m_C)
    {
        m_C.set_frame(frame);
    }
    if (m_O)
    {
        m_O.set_frame(frame);
    }
}

bool mol::SSResidue::is_amino_acid() const
{
    return m_CA && m_C && m_O && m_N;
}

mol::SecondaryStructure mol::SSResidue::secondary_structure() const
{
    return m_structure;
}

void mol::SSResidue::set_secondary_structure(SecondaryStructure const structure)
{
    m_structure = structure;
}

bool mol::SSResidue::is_proline() const
{
    return m_is_proline;
}

bool mol::SSResidue::is_chain_break() const
{
    return m_is_chain_break;
}

void mol::SSResidue::set_chain_break(bool const is_break)
{
    m_is_chain_break = is_break;
}

std::string const& mol::SSResidue::chain() const
{
    return m_chain;
}

mol::Atom const& mol::SSResidue::N() const
{
    return m_N;
}

mol::Atom const& mol::SSResidue::CA() const
{
    return m_CA;
}

mol::Atom const& mol::SSResidue::C() const
{
    return m_C;
}

mol::Atom const& mol::SSResidue::O() const
{
    return m_O;
}

mol::DSSP::DSSP(MolSystem const& molecule)
: m_num_prot_aa{0}
{
    ResidueSel residues(molecule.atoms(0));
    for (mol::Residue residue : residues)
    {
        SSResidue& ss_res = m_residues.emplace_back(residue);
        if (ss_res.is_amino_acid())
        {
            m_num_prot_aa++;
        }

    }
}

void mol::DSSP::run(mol::Frame frame)
{
    dssp::Protein protein(m_num_prot_aa);
    for (size_t index = 0; index < m_residues.size(); index++)
    {
        SSResidue& residue = m_residues[index];
        residue.set_frame(frame);
        if (residue.is_amino_acid())
        {
            protein.emplace_residue(index, residue);
        }
    }

    protein.compute_secondary_structure();

    for (dssp::Residue const& residue : protein.residues())
    {
        SSResidue& ss_residue = m_residues[residue.external_index()];
        ss_residue.set_secondary_structure(residue.structure);
        ss_residue.set_chain_break(residue.is_chain_break());
    }
}

std::vector<mol::SSResidue> const& mol::DSSP::residues() const
{
    return m_residues;
}

std::vector<mol::SecondaryStructure> mol::DSSP::secondary_structures() const
{
    std::vector<mol::SecondaryStructure> structures;
    structures.reserve(m_residues.size());

    for (SSResidue const& residue : m_residues)
    {
            structures.push_back(residue.secondary_structure());
    }

    return structures;
}
