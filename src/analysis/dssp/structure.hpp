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
    MProtein(size_t const num_residues)
    {
        m_residues.reserve(num_residues);
    }
    void CalculateSecondaryStructure(bool inPreferPiHelices = true);
    template<class... Args>
    MResidue& emplace_residue(Args&&... args)
    {
        return m_residues.emplace_back(m_residues.size(), std::forward<Args>(args)...);
    }
    bool is_last(size_t const index) const {return index == m_residues.size() - 1;}

private:
    // TODO revise all arguments to use references whenever possible
    void CalculateHBondEnergies();
    void CalculateAlphaHelices(bool inPreferPiHelices);
    void CalculateBetaSheets();
    bool test_bond(size_t const first, size_t const second) const;
    MBridgeType test_bridge(size_t const first, size_t const second) const;
    double calculate_Hbond_energy(MResidue& inDonor, MResidue& inAcceptor);
    bool no_chain_break(size_t const fisrt, size_t const last) const;
    double compute_kappa(size_t const index) const;

    std::vector<MResidue> m_residues;
};

struct HBond {
    HBond()
    : residue{0}
    , energy{0.0}
    {}

    size_t residue;
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
    MResidue() = default;
    MResidue(size_t const index, std::string chain_id, bool const isProline, MResidue* previous, mol::Point3 N, mol::Point3 CA, mol::Point3 C, mol::Point3 O);
    MResidue(MResidue&& other) = default;
    MResidue(MResidue const&) = delete;
    MResidue& operator=(MResidue const&) = delete;

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

    void compute_H();
    MHelixFlag helix_flag(uint32_t inHelixStride) const;
    void set_helix_flag(uint32_t inHelixStride, MHelixFlag inHelixFlag);
    bool is_helix_start(uint32_t inHelixStride) const;
    bool is_valid_distance(MResidue const* inNext) const;

    size_t index;
    std::string chain_id;
    MResidue* previous;
    MResidue* next;
    mol::SecondaryStructure structure;
    HBond h_bond_donor[2], h_bond_acceptor[2];
    MBridgePartner beta_partner[2];
    uint32_t sheet;
    bool is_bend;
    bool is_proline;
    bool is_chain_break;

private: // TODO make everything public? And a struct?
    mol::Point3 m_N, m_CA, m_C, m_O, m_H;
    MHelixFlag m_helix_flags[3];
};

double distance(mol::Point3 const& first, mol::Point3 const& second);

} // namespace dssp
