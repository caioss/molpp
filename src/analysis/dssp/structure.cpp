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

#include "analysis/dssp/structure.hpp"
#include "analysis/dssp/DSSP.hpp"
#include <cstdint>
#include <iostream>
#include <set>
#include <functional>
#include <list>
#include <cassert>
#include "structure.hpp"

double constexpr kPI = 4 * std::atan(1.0),
                 kSSBridgeDistance = 3.0,
                 kMinimalDistance = 0.5,
                 kMinimalCADistance = 9.0,
                 kMinHBondEnergy = -9.9,
                 kMaxHBondEnergy = -0.5,
                 kCouplingConstant = -332 * 0.42 * 0.2,
                 kMaxPeptideBondLength = 2.5;

double CosinusAngle(const mol::Point3& atom1, const mol::Point3& atom2, const mol::Point3& atom3, const mol::Point3& atom4)
{
    mol::Point3 v12 = atom1 - atom2;
    mol::Point3 v34 = atom3 - atom4;

    double const x = v12.dot(v12) * v34.dot(v34);
    if (x > 0)
    {
        return v12.dot(v34) / sqrt(x);
    }

    return 0;
}

double dssp::distance(mol::Point3 const& first, mol::Point3 const& second)
{
    return (first - second).norm();
}

bool dssp::MProtein::test_bond(dssp::MResidue const* first, dssp::MResidue const* second)
{
    return (first->h_bond_acceptor[0].residue == second && first->h_bond_acceptor[0].energy < kMaxHBondEnergy) || (first->h_bond_acceptor[1].residue == second && first->h_bond_acceptor[1].energy < kMaxHBondEnergy);
}

dssp::MBridgeType dssp::MProtein::test_bridge(dssp::MResidue const* first, dssp::MResidue const* second)
{ // I. a d II. a d parallel
    dssp::MResidue const* a = first->previous; // \ /
    dssp::MResidue const* b = first; // b e b e
    dssp::MResidue const* c = first->next; // / \ ..
    dssp::MResidue const* d = second->previous; // c f c f
    dssp::MResidue const* e = second; //
    dssp::MResidue const* f = second->next; // III. a <- f IV. a f antiparallel
                                                      //
    dssp::MBridgeType result = btNoBridge; // b e b <-> e
                                     //
                                     // c -> d c d

    if (a && c && no_chain_break(a, c) && d && f && no_chain_break(d, f))
    {
        if ((test_bond(c, e) && test_bond(e, a)) || (test_bond(f, b) && test_bond(b, d)))
            result = btParallel;
        else if ((test_bond(c, d) && test_bond(f, a)) || (test_bond(e, b) && test_bond(b, e)))
            result = btAntiParallel;
    }

    return result;
}

double dssp::MProtein::calculate_Hbond_energy(dssp::MResidue* inDonor, dssp::MResidue* inAcceptor)
{
    double result = 0;

    if (!inDonor->is_proline)
    {
        double distanceHO = dssp::distance(inDonor->H(), inAcceptor->O());
        double distanceHC = dssp::distance(inDonor->H(), inAcceptor->C());
        double distanceNC = dssp::distance(inDonor->N(), inAcceptor->C());
        double distanceNO = dssp::distance(inDonor->N(), inAcceptor->O());

        if (distanceHO < kMinimalDistance || distanceHC < kMinimalDistance || distanceNC < kMinimalDistance || distanceNO < kMinimalDistance)
        {
            result = kMinHBondEnergy;
        }
        else
        {
            result = kCouplingConstant / distanceHO - kCouplingConstant / distanceHC + kCouplingConstant / distanceNC - kCouplingConstant / distanceNO;
        }

        if (result < kMinHBondEnergy)
        {
            result = kMinHBondEnergy;
        }
    }

    // update donor
    if (result < inDonor->h_bond_acceptor[0].energy)
    {
        inDonor->h_bond_acceptor[1] = inDonor->h_bond_acceptor[0];
        inDonor->h_bond_acceptor[0].residue = inAcceptor;
        inDonor->h_bond_acceptor[0].energy = result;
    }
    else if (result < inDonor->h_bond_acceptor[1].energy)
    {
        inDonor->h_bond_acceptor[1].residue = inAcceptor;
        inDonor->h_bond_acceptor[1].energy = result;
    }

    // and acceptor
    if (result < inAcceptor->h_bond_donor[0].energy)
    {
        inAcceptor->h_bond_donor[1] = inAcceptor->h_bond_donor[0];
        inAcceptor->h_bond_donor[0].residue = inDonor;
        inAcceptor->h_bond_donor[0].energy = result;
    }
    else if (result < inAcceptor->h_bond_donor[1].energy)
    {
        inAcceptor->h_bond_donor[1].residue = inDonor;
        inAcceptor->h_bond_donor[1].energy = result;
    }

    return result;
}

bool dssp::MProtein::no_chain_break(dssp::MResidue const* from, dssp::MResidue const* to) const
{
    for (dssp::MResidue const* r = from; r != to; r = r->next)
    {
        dssp::MResidue const* next = r->next;
        if (!next || next->is_chain_break)
        {
            return false;
        }
    }
    return true;
}

double dssp::MProtein::compute_kappa(MResidue const* residue) const
{
    const dssp::MResidue* prevPrev = residue->previous ? residue->previous->previous : nullptr;
    const dssp::MResidue* nextNext = residue->next ? residue->next->next : nullptr;

    if (prevPrev != nullptr && nextNext != nullptr && no_chain_break(prevPrev, nextNext))
    {
        double ckap = CosinusAngle(residue->CA(), prevPrev->CA(), nextNext->CA(), residue->CA());
        double skap = sqrt(1 - ckap * ckap);
        return atan2(skap, ckap) * 180 / kPI;
    }
    else
    {
        return 360;
    }
}

struct MBridge
{
    dssp::MBridgeType type;
    uint32_t sheet, ladder;
    std::set<MBridge*> link;
    std::list<size_t> i, j;
    std::string chainI; // TODO non-reference is bad
    std::string chainJ; // TODO non-reference is bad

    MBridge(dssp::MBridgeType type, std::string const& I, std::string const& J)
    : type{type}
    , chainI{I}
    , chainJ{J}
    {
    }

    bool operator<(const MBridge& b) const
    {
        return chainI < b.chainI || (chainI == b.chainI && i.front() < b.i.front());
    }
};

dssp::MResidue::MResidue(std::string chain_id, bool const isProline, MResidue* previous, mol::Point3 N, mol::Point3 CA, mol::Point3 C, mol::Point3 O)
: chain_id(chain_id)
, previous{previous}
, next{nullptr}
, structure(mol::SecondaryStructure::Loop)
, sheet(0)
, is_proline(isProline)
, is_chain_break{false}
, m_N{N}
, m_CA{CA}
, m_C{C}
, m_O{O}
{
    std::fill(m_helix_flags, m_helix_flags + 3, helixNone);
    compute_H();

    if (previous)
    {
        previous->next = this;
    }

    if (previous && !previous->is_valid_distance(this))
    {
        is_chain_break = true;
    }
}

void dssp::MResidue::compute_H()
{
    m_H = m_N;

    if (!is_proline && previous)
    {
        m_H += (previous->C() - previous->O()).normalized();
    }
}

bool dssp::MResidue::is_valid_distance(dssp::MResidue const* inNext) const // TODO accept a reference
{
    return distance(C(), inNext->N()) <= kMaxPeptideBondLength;
}

dssp::MHelixFlag dssp::MResidue::helix_flag(uint32_t inHelixStride) const
{
    assert(inHelixStride == 3 || inHelixStride == 4 || inHelixStride == 5);
    return m_helix_flags[inHelixStride - 3];
}

bool dssp::MResidue::is_helix_start(uint32_t inHelixStride) const
{
    assert(inHelixStride == 3 || inHelixStride == 4 || inHelixStride == 5);
    return m_helix_flags[inHelixStride - 3] == helixStart || m_helix_flags[inHelixStride - 3] == helixStartAndEnd;
}

void dssp::MResidue::set_helix_flag(uint32_t inHelixStride, dssp::MHelixFlag inHelixFlag)
{
    assert(inHelixStride == 3 || inHelixStride == 4 || inHelixStride == 5);
    m_helix_flags[inHelixStride - 3] = inHelixFlag;
}

void dssp::MProtein::CalculateSecondaryStructure(bool inPreferPiHelices)
{
    CalculateHBondEnergies();
    CalculateBetaSheets();
    CalculateAlphaHelices(inPreferPiHelices);
}

void dssp::MProtein::CalculateHBondEnergies()
{
    // Calculate the HBond energies
    for (uint32_t i = 0; i + 1 < m_residues.size(); ++i)
    {
        dssp::MResidue* ri = m_residues[i];

        for (uint32_t j = i + 1; j < m_residues.size(); ++j)
        {
            dssp::MResidue* rj = m_residues[j];

            if (distance(ri->CA(), rj->CA()) < kMinimalCADistance)
            {
                calculate_Hbond_energy(ri, rj);
                if (j != i + 1)
                {
                    calculate_Hbond_energy(rj, ri);
                }
            }
        }
    }
}

void dssp::MProtein::CalculateAlphaHelices(bool inPreferPiHelices)
{
    std::vector<dssp::MResidue*>& res = m_residues;
    for (uint32_t stride = 3; stride <= 5; ++stride)
    {
        if (res.size() < stride)
        {
            continue;
        }

        for (uint32_t i = 0; i + stride < res.size(); ++i)
        {
            if (test_bond(res[i + stride], res[i]) && no_chain_break(res[i], res[i + stride]))
            {
                res[i + stride]->set_helix_flag(stride, helixEnd);
                for (uint32_t j = i + 1; j < i + stride; ++j)
                {
                    if (res[j]->helix_flag(stride) == helixNone)
                        res[j]->set_helix_flag(stride, helixMiddle);
                }

                if (res[i]->helix_flag(stride) == helixEnd)
                {
                    res[i]->set_helix_flag(stride, helixStartAndEnd);
                }
                else
                {
                    res[i]->set_helix_flag(stride, helixStart);
                }
            }
        }
    }

    for (dssp::MResidue* r : m_residues)
    {
        double kappa = compute_kappa(r);
        r->is_bend = kappa != 360 && kappa > 70;
    }

    for (uint32_t i = 1; i + 4 < m_residues.size(); ++i)
    {
        if (m_residues[i]->is_helix_start(4) && m_residues[i - 1]->is_helix_start(4))
        {
            for (uint32_t j = i; j <= i + 3; ++j)
            {
                m_residues[j]->structure = mol::SecondaryStructure::Helix;
            }
        }
    }

    for (uint32_t i = 1; i + 3 < m_residues.size(); ++i)
    {
        if (m_residues[i]->is_helix_start(3) && m_residues[i - 1]->is_helix_start(3))
        {
            bool empty = true;
            for (uint32_t j = i; empty && j <= i + 2; ++j)
            {
                empty = m_residues[j]->structure == mol::SecondaryStructure::Loop || m_residues[j]->structure == mol::SecondaryStructure::Helix3;
            }

            if (empty)
            {
                for (uint32_t j = i; j <= i + 2; ++j)
                {
                    m_residues[j]->structure = mol::SecondaryStructure::Helix3;
                }
            }
        }
    }

    for (uint32_t i = 1; i + 5 < m_residues.size(); ++i)
    {
        if (m_residues[i]->is_helix_start(5) && m_residues[i - 1]->is_helix_start(5))
        {
            bool empty = true;
            for (uint32_t j = i; empty && j <= i + 4; ++j)
                empty = m_residues[j]->structure == mol::SecondaryStructure::Loop || m_residues[j]->structure == mol::SecondaryStructure::Helix5 ||
                        (inPreferPiHelices && m_residues[j]->structure == mol::SecondaryStructure::Helix);
            if (empty)
            {
                for (uint32_t j = i; j <= i + 4; ++j)
                    m_residues[j]->structure = mol::SecondaryStructure::Helix5;
            }
        }
    }

    for (uint32_t i = 1; i + 1 < m_residues.size(); ++i)
    {
        if (m_residues[i]->structure == mol::SecondaryStructure::Loop)
        {
            bool isTurn = false;
            for (uint32_t stride = 3; stride <= 5 && not isTurn; ++stride)
            {
                for (uint32_t k = 1; k < stride && not isTurn; ++k)
                    isTurn = (i >= k) && m_residues[i - k]->is_helix_start(stride);
            }

            if (isTurn)
                m_residues[i]->structure = mol::SecondaryStructure::Turn;
            else if (m_residues[i]->is_bend)
                m_residues[i]->structure = mol::SecondaryStructure::Bend;
        }
    }
}

void dssp::MProtein::CalculateBetaSheets()
{
    // Calculate Bridges
    std::vector<MBridge> bridges;
    if (m_residues.size() > 4)
    {
        for (size_t i = 1; i + 4 < m_residues.size(); ++i)
        {
            dssp::MResidue* ri = m_residues[i];

            for (size_t j = i + 3; j + 1 < m_residues.size(); ++j)
            {
                dssp::MResidue* rj = m_residues[j];

                dssp::MBridgeType type = test_bridge(ri, rj);
                if (type == btNoBridge)
                    continue;

                bool found = false;
                for (MBridge& bridge : bridges)
                {
                    if (type != bridge.type || i != bridge.i.back() + 1)
                        continue;

                    if (type == btParallel && bridge.j.back() + 1 == j)
                    {
                        bridge.i.push_back(i);
                        bridge.j.push_back(j);
                        found = true;
                        break;
                    }

                    if (type == btAntiParallel && bridge.j.front() - 1 == j)
                    {
                        bridge.i.push_back(i);
                        bridge.j.push_front(j);
                        found = true;
                        break;
                    }
                }

                if (!found)
                {
                    MBridge bridge(type, ri->chain_id, rj->chain_id);
                    bridge.i.push_back(i);
                    bridge.j.push_back(j);
                    bridges.push_back(bridge);
                }
            }
        }
    }

    // extend ladders
    std::sort(bridges.begin(), bridges.end());

    for (uint32_t i = 0; i < bridges.size(); ++i)
    {
        for (uint32_t j = i + 1; j < bridges.size(); ++j)
        {
            uint32_t ibi = bridges[i].i.front();
            uint32_t iei = bridges[i].i.back();
            uint32_t jbi = bridges[i].j.front();
            uint32_t jei = bridges[i].j.back();
            uint32_t ibj = bridges[j].i.front();
            uint32_t iej = bridges[j].i.back();
            uint32_t jbj = bridges[j].j.front();
            uint32_t jej = bridges[j].j.back();

            if (bridges[i].type != bridges[j].type ||
                no_chain_break(m_residues[std::min(ibi, ibj)], m_residues[std::max(iei, iej)]) == false ||
                no_chain_break(m_residues[std::min(jbi, jbj)], m_residues[std::max(jei, jej)]) == false ||
                ibj - iei >= 6 ||
                (iei >= ibj && ibi <= iej))
            {
                continue;
            }

            bool bulge;
            if (bridges[i].type == btParallel)
            {
                bulge = ((jbj - jei < 6 && ibj - iei < 3) || (jbj - jei < 3));
            }
            else
            {
                bulge = ((jbi - jej < 6 && ibj - iei < 3) || (jbi - jej < 3));
            }

            if (bulge)
            {
                bridges[i].i.insert(bridges[i].i.end(), bridges[j].i.begin(), bridges[j].i.end());
                if (bridges[i].type == btParallel)
                {
                    bridges[i].j.insert(bridges[i].j.end(), bridges[j].j.begin(), bridges[j].j.end());
                }
                else
                {
                    bridges[i].j.insert(bridges[i].j.begin(), bridges[j].j.begin(), bridges[j].j.end());
                }

                bridges.erase(bridges.begin() + j);
                --j;
            }
        }
    }

    for (MBridge& bridge : bridges)
    {
        // find out if any of the i and j set members already have
        // a bridge assigned, if so, we're assigning bridge 2
        mol::SecondaryStructure ss = mol::SecondaryStructure::Bridge;
        if (bridge.i.size() > 1)
            ss = mol::SecondaryStructure::Strand;

        for (size_t i = bridge.i.front(); i <= bridge.i.back(); ++i)
        {
            if (m_residues[i]->structure != mol::SecondaryStructure::Strand)
            {
                m_residues[i]->structure = ss;
            }
            m_residues[i]->sheet = bridge.sheet;
        }

        for (size_t i = bridge.j.front(); i <= bridge.j.back(); ++i)
        {
            if (m_residues[i]->structure != mol::SecondaryStructure::Strand)
            {
                m_residues[i]->structure = ss;
            }
            m_residues[i]->sheet = bridge.sheet;
        }
    }
}
