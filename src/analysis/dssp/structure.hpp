/*
 * Copyright Maarten L. Hekkelman, Radboud University 2008-2011.
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *
 * This is a port of version 2.2.1 of DSSP to Molpp. DSSP is
 * distributed under the Boost Software License, Version 1.0.
 * Modifications here present were made by Caio S. Souza and are
 * distributed under the BSD 3-Clause License.
 */

// TODO
#pragma once

// TODO rename and separate classes

#include <molpp/MolppCore.hpp>
#include <molpp/Atom.hpp> // TODO needed
#include "analysis/dssp/DSSP.hpp"
#include <string>
#include <vector>

namespace dssp
{

enum MBridgeType {
    btNoBridge, btParallel, btAntiParallel
};

enum MHelixFlag {
    helixNone, helixStart, helixEnd, helixStartAndEnd, helixMiddle
};

class MResidue;
class MChain;

// TODO merge with the analysis DSSP?
class MProtein
{
public:
    void CalculateSecondaryStructure(bool inPreferPiHelices = true);
    void add_residue(MResidue* residue) {m_residues.push_back(residue);}
    std::vector<MResidue*> const& residues() {return m_residues;}

private:
    void CalculateHBondEnergies();
    void CalculateAlphaHelices(bool inPreferPiHelices);
    void CalculateBetaSheets();
    bool test_bond(MResidue const* first, MResidue const* second);
    MBridgeType test_bridge(MResidue const* first, MResidue const* second);
    double calculate_Hbond_energy(MResidue* inDonor, MResidue* inAcceptor);
    bool no_chain_break(MResidue const* from, MResidue const* to) const;
    double compute_kappa(MResidue const* residue) const;

    std::vector<MResidue*> m_residues;
};

struct HBond {
    HBond()
    : residue{nullptr}
    , energy{0.0}
    {}

    MResidue const* residue;
    double energy;
};

struct MBridgePartner {
    MBridgePartner()
    : residue{nullptr}
    {}

    MResidue const* residue;
    uint32_t ladder;
    bool parallel;
};

class MResidue
{
public:
    MResidue(std::string chain_id, bool const isProline);

    mol::Point3 const& N() const
    {
        return m_N;
    }

    mol::Point3 const& CA() const
    {
        return m_CA;
    }

    mol::Point3 const& C() const
    {
        return m_C;
    }

    mol::Point3 const& O() const
    {
        return m_O;
    }

    mol::Point3 const& H() const
    {
        return m_H;
    }

    void update_positions(mol::Atom const& N, mol::Atom const& CA, mol::Atom const& C, mol::Atom const& O);
    MHelixFlag helix_flag(uint32_t inHelixStride) const;
    void set_helix_flag(uint32_t inHelixStride, MHelixFlag inHelixFlag);
    bool is_helix_start(uint32_t inHelixStride) const;
    bool is_valid_distance(MResidue const* inNext) const;

    std::string chain_id;
    MResidue const* previous;
    MResidue const* next;
    mol::SecondaryStructure structure;
    HBond h_bond_donor[2], h_bond_acceptor[2];
    MBridgePartner beta_partner[2];
    uint32_t sheet;
    bool is_bend;
    bool is_proline;
    bool is_chain_break;

private:
    mol::Point3 m_N, m_CA, m_C, m_O;
    mol::Point3 m_H;
    MHelixFlag m_helix_flags[3];
};

double distance(mol::Point3 const& first, mol::Point3 const& second);

} // namespace dssp
