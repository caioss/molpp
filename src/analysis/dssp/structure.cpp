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
#include <set>
#include <list>

double constexpr PI = 4 * std::atan(1.0);
double constexpr SS_BRIDGE_DISTANCE = 3.0;
double constexpr MIN_DISTANCE = 0.5;
double constexpr MIN_CA_DISTANCE = 9.0;
double constexpr MIN_HBOND_ENERGY = -9.9;
double constexpr MAX_HBOND_ENERGY = -0.5;
double constexpr COUPLING_CONSTANT = -332 * 0.42 * 0.2;
double constexpr MAX_PEPTIDE_BOND_LENGTH = 2.5;

double cosinus_angle(const mol::Point3& atom1, const mol::Point3& atom2, const mol::Point3& atom3, const mol::Point3& atom4)
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

bool dssp::MProtein::test_bond(size_t const first, size_t const second) const
{
    MResidue const& first_residue = m_residues[first];
    HBond const& hb1 = first_residue.h_bond_acceptor[0];
    HBond const& hb2 = first_residue.h_bond_acceptor[1];

    return (hb1.residue == second && hb1.energy < MAX_HBOND_ENERGY) || (hb2.residue == second && hb2.energy < MAX_HBOND_ENERGY);
}

dssp::BridgeType dssp::MProtein::test_bridge(size_t const first, size_t const second) const
{
    if (first == 0 || is_last(first) || second == 0 || is_last(second))
    {
        return BridgeType::None;
    }

    // I. a d II. a d parallel
    //    \ /
    //    b e b e
    //    / \ ..
    //    c f c f
    //
    // III. a <- f IV. a f antiparallel
    //
    // b e b <-> e
    //
    // c -> d c d

    size_t const a = first - 1;
    size_t const b = first;
    size_t const c = first + 1;
    size_t const d = second - 1;
    size_t const e = second;
    size_t const f = second + 1;

    if (no_chain_break(a, c) && no_chain_break(d, f))
    {
        if ((test_bond(c, e) && test_bond(e, a)) || (test_bond(f, b) && test_bond(b, d)))
        {
            return BridgeType::Parallel;
        }
        else if ((test_bond(c, d) && test_bond(f, a)) || (test_bond(e, b) && test_bond(b, e)))
        {
            return BridgeType::AntiParallel;
        }
    }

    return BridgeType::None;
}

double dssp::MProtein::compute_h_bond(dssp::MResidue& donor, dssp::MResidue& acceptor)
{
    double result = 0;

    if (!donor.is_proline)
    {
        double const dist_H_O = dssp::distance(donor.H(), acceptor.O());
        double const dist_H_C = dssp::distance(donor.H(), acceptor.C());
        double const dist_N_C = dssp::distance(donor.N(), acceptor.C());
        double const dist_N_O = dssp::distance(donor.N(), acceptor.O());

        if (dist_H_O < MIN_DISTANCE || dist_H_C < MIN_DISTANCE || dist_N_C < MIN_DISTANCE || dist_N_O < MIN_DISTANCE)
        {
            result = MIN_HBOND_ENERGY;
        }
        else
        {
            result = COUPLING_CONSTANT / dist_H_O - COUPLING_CONSTANT / dist_H_C + COUPLING_CONSTANT / dist_N_C - COUPLING_CONSTANT / dist_N_O;
        }

        if (result < MIN_HBOND_ENERGY)
        {
            result = MIN_HBOND_ENERGY;
        }
    }

    // Update donor
    if (result < donor.h_bond_acceptor[0].energy)
    {
        donor.h_bond_acceptor[1] = donor.h_bond_acceptor[0];
        donor.h_bond_acceptor[0].residue = acceptor.index;
        donor.h_bond_acceptor[0].energy = result;
    }
    else if (result < donor.h_bond_acceptor[1].energy)
    {
        donor.h_bond_acceptor[1].residue = acceptor.index;
        donor.h_bond_acceptor[1].energy = result;
    }

    // Update acceptor
    if (result < acceptor.h_bond_donor[0].energy)
    {
        acceptor.h_bond_donor[1] = acceptor.h_bond_donor[0];
        acceptor.h_bond_donor[0].residue = donor.index;
        acceptor.h_bond_donor[0].energy = result;
    }
    else if (result < acceptor.h_bond_donor[1].energy)
    {
        acceptor.h_bond_donor[1].residue = donor.index;
        acceptor.h_bond_donor[1].energy = result;
    }

    return result;
}

bool dssp::MProtein::no_chain_break(size_t const fisrt, size_t const last) const
{
    for (size_t index = fisrt; index < last; index++)
    {
        if (is_last(index) || m_residues[index + 1].is_chain_break)
        {
            return false;
        }
    }
    return true;
}

double dssp::MProtein::compute_kappa(size_t const index) const
{
    if (index < 2 || index > m_residues.size() - 2)
    {
        return 360;
    }

    size_t const prev_prev = index - 2;
    size_t const next_next = index + 2;

    if (!no_chain_break(prev_prev, next_next))
    {
        return 360;
    }

    double const ckap = cosinus_angle(m_residues[index].CA(), m_residues[prev_prev].CA(), m_residues[next_next].CA(), m_residues[index].CA());
    double const skap = sqrt(1 - ckap * ckap);
    return atan2(skap, ckap) * 180 / PI;
}

bool dssp::MProtein::is_last(size_t const index) const
{
    return index == m_residues.size() - 1;
}

struct MBridge
{
    dssp::BridgeType type;
    uint32_t sheet, ladder;
    std::set<MBridge*> link;
    std::list<size_t> i, j;
    std::string chainI; // TODO non-reference is bad
    std::string chainJ; // TODO non-reference is bad

    MBridge(dssp::BridgeType type, std::string const& I, std::string const& J)
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

dssp::MResidue::MResidue(size_t const index, std::string chain_id, bool const is_proline, mol::Point3 N, mol::Point3 CA, mol::Point3 C, mol::Point3 O)
: index{index}
, chain_id(chain_id)
, structure(mol::Loop)
, sheet(0)
, is_proline(is_proline)
, is_chain_break{false}
, m_N{N}
, m_CA{CA}
, m_C{C}
, m_O{O}
, m_H{N}
{
    std::fill(m_helix_flags, m_helix_flags + 3, HelixType::None);
}

void dssp::MResidue::set_previous(MResidue const& previous)
{
    if (is_proline)
    {
        return;
    }

    m_H += (previous.C() - previous.O()).normalized();

    if (!is_valid_distance(previous))
    {
        is_chain_break = true;
    }
}

bool dssp::MResidue::is_valid_distance(dssp::MResidue const& next) const
{
    return distance(N(), next.C()) <= MAX_PEPTIDE_BOND_LENGTH;
}

dssp::HelixType dssp::MResidue::helix_flag(uint32_t const stride) const
{
    return m_helix_flags[stride - 3];
}

bool dssp::MResidue::is_helix_start(uint32_t const stride) const
{
    return m_helix_flags[stride - 3] == HelixType::Start || m_helix_flags[stride - 3] == HelixType::StartAndEnd;
}

void dssp::MResidue::set_helix_flag(uint32_t const stride, dssp::HelixType const type)
{
    m_helix_flags[stride - 3] = type;
}

dssp::MProtein::MProtein(size_t const num_residues)
{
    m_residues.reserve(num_residues);
}

void dssp::MProtein::compute_secondary_structure(bool prefer_pi_helices)
{
    compute_h_bond_energies();
    compute_sheets();
    compute_helices(prefer_pi_helices);
}

void dssp::MProtein::compute_h_bond_energies()
{
    size_t const num_residues = m_residues.size();
    for (uint32_t i = 0; i + 1 < num_residues; ++i)
    {
        dssp::MResidue& ri = m_residues[i];

        for (uint32_t j = i + 1; j < num_residues; ++j)
        {
            dssp::MResidue& rj = m_residues[j];

            if (distance(ri.CA(), rj.CA()) < MIN_CA_DISTANCE)
            {
                compute_h_bond(ri, rj);
                if (j != i + 1)
                {
                    compute_h_bond(rj, ri);
                }
            }
        }
    }
}

void dssp::MProtein::compute_helices(bool prefer_pi_helices)
{
    for (size_t stride = 3; stride <= 5; ++stride)
    {
        if (m_residues.size() < stride)
        {
            continue;
        }

        for (size_t i = 0; i + stride < m_residues.size(); ++i)
        {
            MResidue& current = m_residues[i];

            if (test_bond(i + stride, i) && no_chain_break(i, i + stride))
            {
                m_residues[i + stride].set_helix_flag(stride, HelixType::End);
                for (size_t j = i + 1; j < i + stride; ++j)
                {
                    MResidue& between = m_residues[j];
                    if (between.helix_flag(stride) == HelixType::None)
                    {
                        between.set_helix_flag(stride, HelixType::Middle);
                    }
                }

                if (current.helix_flag(stride) == HelixType::End)
                {
                    current.set_helix_flag(stride, HelixType::StartAndEnd);
                }
                else
                {
                    current.set_helix_flag(stride, HelixType::Start);
                }
            }
        }
    }

    for (dssp::MResidue& residue : m_residues)
    {
        double kappa = compute_kappa(residue.index);
        residue.is_bend = kappa != 360 && kappa > 70;
    }

    size_t const num_residues = m_residues.size();
    for (size_t i = 1; i < num_residues; ++i)
    {
        if (m_residues[i].is_helix_start(4) && m_residues[i - 1].is_helix_start(4))
        {
            for (size_t j = i; j <= i + 3; ++j)
            {
                m_residues[j].structure = mol::Helix;
            }
        }
    }

    for (size_t i = 1; i < num_residues; ++i)
    {
        if (m_residues[i].is_helix_start(3) && m_residues[i - 1].is_helix_start(3))
        {
            bool empty = true;
            for (size_t j = i; empty && j <= i + 2; ++j)
            {
                mol::SecondaryStructure const structure = m_residues[j].structure;
                empty = structure == mol::Loop || structure == mol::Helix3;
            }

            if (empty)
            {
                for (size_t j = i; j <= i + 2; ++j)
                {
                    m_residues[j].structure = mol::Helix3;
                }
            }
        }
    }

    for (size_t i = 1; i < num_residues; ++i)
    {
        if (m_residues[i].is_helix_start(5) && m_residues[i - 1].is_helix_start(5))
        {
            bool empty = true;
            for (size_t j = i; empty && j <= i + 4; ++j)
            {
                mol::SecondaryStructure const structure = m_residues[j].structure;
                empty = structure == mol::Loop
                     || structure == mol::Helix5
                     || (prefer_pi_helices && structure == mol::Helix);
            }

            if (empty)
            {
                for (size_t j = i; j <= i + 4; ++j)
                {
                    m_residues[j].structure = mol::Helix5;
                }
            }
        }
    }

    for (size_t i = 1; i < num_residues - 1; ++i)
    {
        MResidue& residue = m_residues[i];
        if (residue.structure == mol::Loop)
        {
            bool is_turn = false;
            for (size_t stride = 3; stride <= 5 && !is_turn; ++stride)
            {
                for (size_t k = 1; k < stride && !is_turn; ++k)
                {
                    is_turn = (i >= k) && m_residues[i - k].is_helix_start(stride);
                }
            }

            if (is_turn)
            {
                residue.structure = mol::Turn;
            }
            else if (residue.is_bend)
            {
                residue.structure = mol::Bend;
            }
        }
    }
}

void dssp::MProtein::compute_sheets()
{
    // Calculate Bridges
    std::vector<MBridge> bridges;
    if (m_residues.size() > 4)
    {
        for (size_t i = 1; i + 4 < m_residues.size(); ++i)
        {
            dssp::MResidue& ri = m_residues[i];

            for (size_t j = i + 3; j + 1 < m_residues.size(); ++j)
            {
                dssp::MResidue& rj = m_residues[j];

                dssp::BridgeType type = test_bridge(i, j);
                if (type == BridgeType::None)
                {
                    continue;
                }

                bool found = false;
                for (MBridge& bridge : bridges)
                {
                    if (type != bridge.type || i != bridge.i.back() + 1)
                    {
                        continue;
                    }

                    if (type == BridgeType::Parallel && bridge.j.back() + 1 == j)
                    {
                        bridge.i.push_back(i);
                        bridge.j.push_back(j);
                        found = true;
                        break;
                    }

                    if (type == BridgeType::AntiParallel && bridge.j.front() - 1 == j)
                    {
                        bridge.i.push_back(i);
                        bridge.j.push_front(j);
                        found = true;
                        break;
                    }
                }

                if (!found)
                {
                    MBridge bridge(type, ri.chain_id, rj.chain_id);
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

            if (bridges[i].type != bridges[j].type
                || no_chain_break(std::min(ibi, ibj), std::max(iei, iej)) == false
                || no_chain_break(std::min(jbi, jbj), std::max(jei, jej)) == false
                || ibj - iei >= 6
                || (iei >= ibj && ibi <= iej))
            {
                continue;
            }

            bool bulge;
            if (bridges[i].type == BridgeType::Parallel)
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
                if (bridges[i].type == BridgeType::Parallel)
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
        mol::SecondaryStructure ss = mol::Bridge;
        if (bridge.i.size() > 1)
        {
            ss = mol::Strand;
        }

        for (size_t i = bridge.i.front(); i <= bridge.i.back(); ++i)
        {
            MResidue& residue = m_residues[i];
            if (residue.structure != mol::Strand)
            {
                residue.structure = ss;
            }
            residue.sheet = bridge.sheet;
        }

        for (size_t i = bridge.j.front(); i <= bridge.j.back(); ++i)
        {
            MResidue& residue = m_residues[i];
            if (residue.structure != mol::Strand)
            {
                residue.structure = ss;
            }
            residue.sheet = bridge.sheet;
        }
    }
}
